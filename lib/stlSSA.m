function [Xfinal,R,IntG,StepsTillTend] = stlSSA(X0, P, Tend) %#codegen
    
    [Sims, Species] = size(X0);
    Rs = size(P,2);
    
    T=zeros(Sims,1);
    R=zeros(Sims,Rs);
    IntG=zeros(Sims,Rs);
    StepsTillTend=zeros(Sims,1);
    X = X0;
    [S, Sed] = getStoichMatrix();

    
    active = T<Tend; 
    activeSims = sum(active);
    Xact = X(active,:);
    Pact = P(active,:);
    Tact = T(active);
    Ract = R(active,:);
    IntGact = IntG(active,:);
    props = computePropensities(Pact,Xact);
    intGG = props./Pact; % Valid for mass-action only

    stepCount = 1;
    while any(active)
        tau = exprnd(1./sum(props,2));
        PropCritSums=cumsum(props,2);
        PropCritSums = PropCritSums./repmat(PropCritSums(:,end),[1 Rs]);
        RandCrit = rand([activeSims 1]);
        nextReactions = Rs+1-sum(PropCritSums-repmat(RandCrit,[1 Rs])>0,2);

        % Set/check new time
        tau = min(tau, Tend-Tact);
        Tact = Tact + tau;

        activeCrit = Tact<Tend;
        Xact(activeCrit,:) = Xact(activeCrit,:) + S(:,nextReactions(activeCrit))';
        Rlinidcs = sub2ind(size(Ract), find(activeCrit), nextReactions(activeCrit));
        Ract(Rlinidcs) = Ract(Rlinidcs) + 1;
        
        % Update integral G 
        intGGprev = intGG;
        props = computePropensities(Pact,Xact);
        intGG = props./Pact; % Valid for mass-action only
        
        intGGprevTau = intGGprev.*repmat(tau,[1 Rs]);
        IntGact = IntGact + intGGprevTau;
        
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
        IntGact = IntG(active,:);
        props = props(newFromOldActive,:);
        intGG = intGG(newFromOldActive,:);

        StepsTillTend(T>=Tend & StepsTillTend==0) = stepCount;
        stepCount = stepCount + 1;
    end

    Xfinal = X;
end
