function [Xfinal,R,IntG,StepsTillTend] = stlTauLeaping(X0, P, Tend, fStoichMatrix, fPropensities, ep, Nc) %#codegen
% STLTAULEAPING performs stochastic simulation with explicit, adaptive tau-leaping [1] 
% with X0 (instances-by-species), parameters P (instances-by-parameters) and end time Tend.
%
% fStoichMatrix and fPropensities are function handles: 
% [StoichMatrix,EductMatrix] = fStoichMatrix(),
% [Propensities] = fPropensities(Parameters,SystemState), Input particles-by-params/species
%
% [1] Cao, Y., Gillespie, D. T., & Petzold, L. R. (2007). Adaptive explicit-implicit 
% tau-leaping method with automatic tau selection. Journal of Chemical Physics, 
% 126(2007), 1–27. http://doi.org/10.1063/1.2745299

    if (isempty(Nc))
        Nc = 20;    % Threshold for critical reaction
    end
    
    if (isempty(ep))
        ep = 0.02;   % Error control parameter
    end
    
    [Sims, Species] = size(X0);
    Rs = size(P,2);
    
    T=zeros(Sims,1);
    R=zeros(Sims,Rs);
    IntG=zeros(Sims,Rs);
    StepsTillTend=zeros(Sims,1);
    X = X0;
    [S, Sed] = fStoichMatrix();
    Sed = abs(Sed);
    HOR = max(repmat(sum(Sed, 1),[Species, 1]).*(Sed>0),[],2);
    HEO = max(Sed,[],2);
    
    Sedeff = S; Sedeff(Sedeff>0) = 0; Sedeff = abs(Sedeff);
    active = T<Tend; 
    activeSims = sum(active);
    Xact = X(active,:);
    Pact = P(active,:);
    Tact = T(active);
    Ract = R(active,:);
    Gact = computeG(Xact, HOR, HEO);
    IntGact = IntG(active,:);
    props = fPropensities(Pact,Xact);
    intGG = props./Pact; % Valid for mass-action only

    stepCount = 1;
    while any(active)
        
        %sprintf('%d\n',stepCount)
        Ssim = repmat(permute(S,[3 1 2]),[activeSims,1,1]);
        Ssim2 = repmat(permute(S.^2,[3 1 2]),[activeSims,1,1]);
        Sedeffsim = repmat(permute(Sedeff,[3 1 2]),[activeSims,1,1]);
        
        % Algorithm step (1) in [1], formula (16)
        propsGt0sim = repmat(permute(props~=0,[1 3 2]),[1 Species 1]);
        maxReactSpec=repmat(Xact,[1,1,Rs])./(Sedeffsim.*propsGt0sim);
        maxReactSpec(isnan(maxReactSpec)) = Inf;
        L=reshape(min(maxReactSpec,[], 2),[activeSims, Rs]);
        RCrit = L<Nc;
        RCritOrZero = RCrit | props==0;
        
        % Algorithm step (2) in [1], formulas (9a) (9b)
        propsNonCrit = props; 
        propsNonCrit(RCrit) = 0;
        propsNonCritSpec = repmat(permute(propsNonCrit,[1 3 2]),[1 Species 1]);
        mu = sum(Ssim.*propsNonCritSpec, 3);
        sig = sum(Ssim2.*propsNonCritSpec, 3);
        eductsNonCrit = any(Sedeffsim.*repmat(permute(~RCritOrZero,[1 3 2]), [1 Species 1])>0,3);

        % Algorithm step (2) in [1], formula (8)
        epsXbyG = ep*Xact./Gact;
        epsXbyG(epsXbyG<1) = 1;
        epsXbyG(~eductsNonCrit) = NaN;
        [tau1, idx1] = nanmin([epsXbyG./abs(mu), epsXbyG.^2./sig Tend-Tact],[],2);
        leapEnd = (idx1==Species*2+1);
        
        % Algorithm step (2) in [1], formula (8)
        propsCrit = props.*RCrit;
        tau2 = 1./sum(propsCrit,2).*log(1./rand([activeSims,1]));
        
        [tau, leap] = nanmin([tau1, tau2],[],2);
        leap = (leap==1);
        
        % Algorithm step (5-7) in [1]
        PropCritSums=cumsum(propsCrit(~leap,:),2);
        PropCritSums = PropCritSums./repmat(PropCritSums(:,end),[1 Rs]);
        RandCrit = rand([nnz(~leap) 1]);
        critReactions = Rs+1-sum(PropCritSums-repmat(RandCrit,[1 Rs])>0,2);
        
        XactNConly = Xact;
        if (nnz(~leap)>0)
            Xact(~leap,:) = Xact(~leap,:) + S(:,critReactions)';
            Rlinidcs = sub2ind(size(Ract), find(~leap), critReactions);
            Ract(Rlinidcs) = Ract(Rlinidcs) + 1;
        end
        
        leapLambdas = propsNonCrit.*repmat(tau,[1 Rs]);
        nonCritReactionFired=zeros(activeSims, Rs);
        nonCritIdcs = leapLambdas~=0;
        lambdas = leapLambdas(nonCritIdcs);
        if (isrow(lambdas))
            lambdas = lambdas';
        end
        
        % do a normal approximation to lambdas > 25
        largeLambdasIdcs = lambdas>40;
        largeLambdas = lambdas(largeLambdasIdcs);
        nonCritIdcs2 = find(nonCritIdcs);
        if any(largeLambdasIdcs)
            nonCritReactionFired(nonCritIdcs2(largeLambdasIdcs)) = round(max(normrnd(largeLambdas, sqrt(largeLambdas)), 0));
        end
        if ~all(largeLambdasIdcs)
            nonCritReactionFired(nonCritIdcs2(~largeLambdasIdcs)) = poissrnd(lambdas(~largeLambdasIdcs));
        end
%         nonCritReactionFired(nonCritIdcs) = poissrnd(lambdas); % old version
%         
        nonCritStateChange = sum(repmat(permute(nonCritReactionFired,[1 3 2]),[1 Species 1]).*Ssim,3);
        Xact = Xact + nonCritStateChange;
        XactNConly = XactNConly + nonCritStateChange;
        Ract = Ract + nonCritReactionFired;
        
        % Set new time, set end time exactly
        Tact = Tact + tau;
        Tact(leap & leapEnd) = Tend;

        % Update integral G 
        intGGprev = intGG;
        props = fPropensities(Pact,Xact);
        intGG = props./Pact; % Valid for mass-action only
        
        propsNCchange = fPropensities(Pact,XactNConly);
        intGGNCchange = propsNCchange./Pact; % Valid for mass-action only

        intGGprevTau = intGGprev.*repmat(tau,[1 Rs]);
        intGGavgTau = (intGGprev+intGGNCchange)/2.*repmat(tau,[1 Rs]);
        %IntGact(RCrit) = IntGact(RCrit) + intGGprevTau(RCrit);
        %IntGact(~RCrit) = IntGact(~RCrit) + intGGavgTau(~RCrit);
        IntGact = IntGact + intGGavgTau;
        
        % Copy active instances
        X(active,:) = Xact;
        T(active) = Tact;
        R(active,:) = Ract;
        IntG(active,:) = IntGact;
        
        % Load remaining active ones
        previouslyActive = active;
        active = T<Tend; 
        activeSims = sum(active);
        [~,newFromOldActive] = ismember(find(active), find(previouslyActive));
        
        Xact = X(active,:);
        Pact = P(active,:);
        Tact = T(active);
        Ract = R(active,:);
        Gact = computeG(Xact, HOR, HEO);
        IntGact = IntG(active,:);
        props = props(newFromOldActive,:);
        intGG = intGG(newFromOldActive,:);

        StepsTillTend(T>=Tend & StepsTillTend==0) = stepCount;
        stepCount = stepCount + 1;
    end

    Xfinal = X;
end

function [G] = computeG(X, HOR, HEO)
    G=zeros(size(X));
    HORs = repmat(HOR', [size(X,1), 1]);
    HEOs = repmat(HEO', [size(X,1), 1]);
    G(HORs==1) = 1;
    G(HORs==2 & HEOs < 2) = 2;
    G(HORs==2 & HEOs == 2) = 2+(X(HORs==2 & HEOs == 2)-1).^-1;
    G(HORs==3 & HEOs < 2) = 3;
    G(HORs==3 & HEOs == 2) = 3/2*(2+(X(HORs==3 & HEOs == 2)-1).^-1);
    G(HORs==3 & HEOs == 3) = (3+(X(HORs==3 & HEOs == 3)-1).^-1+2*(X(HORs==3 & HEOs == 3)-2).^-1);
end