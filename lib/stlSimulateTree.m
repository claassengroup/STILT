function [ Cells ] = stlSimulateTree(X0, theta, NGen, times, opts, sbmlModel)
% STLSIMULATETREE simulates a tree of NGen generations.
%
% X0 ... Initial state of initial cell
% theta ... vector of parameters
% NGen ... number of generations
% times ... cell array (cells-by-1) of time-vectors per cell
% opts ... STILT options structure
% sbmlModel ... SBML model
% 
% Returns:
% Cells ... Cells-by-5 cell array. 

    Cs = 2^NGen-1;
    Cells = cell(Cs,4); % X0, T0, X, R, IntG
    Cells{1,1} = X0;
    Cells{1,2} = times{1};
    
    for cIdx = 1:Cs
        [Cells{cIdx,3:5}] = simCell(Cells{cIdx,1}, theta, Cells{cIdx,2}, length(sbmlModel.reaction), opts);
        [ daughter1ID ] =  cIdx*2;
        if (daughter1ID+1 <= Cs)
            [Cells{daughter1ID,1}, Cells{daughter1ID+1,1}] = stlSimulateCellDivision(Cells{cIdx,3}(end,:), opts);
            Cells{daughter1ID,2} = [Cells{cIdx,2}(end);times{daughter1ID}];
            Cells{daughter1ID+1,2} = [Cells{cIdx,2}(end);times{daughter1ID+1}];
        end
    end
  
end


%%
% loop over number of generations
% for n=1:Ngen
%     cells = 2^(n-1):(2^n - 1);
%     % loop over cells in this generation
%     for k=cells
%         if nargin>3
%             [Xout, Pr{k}, r, int_g] = simCell(model, dt, X0{k}, merr, varargin);
%         else
%             [Xout, Pr{k}, r, int_g] = simCell(model, dt, X0{k}, merr);
%         end
%         
%         % setup daughter cells for next generation
%         m1 = binornd(Xout(end, ix_mRNA), 0.5);
%         m2 = Xout(end,ix_mRNA)-m1;
%         p1 = normrnd(Xout(end,ix_protein)*0.5, sqrt(Xout(end,ix_protein))*0.5);
%         p2 = Xout(end,ix_protein)-p1;
%         
%         X0{2*k}     = [Xout(end,ix_DNA), m1, p1];
%         X0{2*k+1}   = [Xout(end,ix_DNA), m2, p2];
%         
%         % times
%         N_k = length(Pr{k});
% 
%         if isempty(times)
%             times{k} = dt*(0:(N_k-1));
%         else
%             times{k} = times{floor(k/2)}(end) + dt*(1:N_k);
%         end
%         % record into the table
%         output = [output; [times{k}', repmat(k, N_k, 1), Pr{k}, Xout ]];        
%     end
%     r_total = r_total+r;
%     int_g_total = int_g_total + int_g;
%     
% end
% 
% if nargout>0
%     varargout{1}=output;
% end
% 
% if nargout>1
%     varargout{2} = r_total;
% end
% 
% if nargout>2
%     varargout{3} = int_g_total;
% end
% 
%     end
% end


function [X, R, IntG] = simCell(X0, theta, times, reactions, opts)
    Ts = length(times);
    Ss = size(X0,2);
    Rs = reactions;
    
    dts = diff(times);
    X = zeros(Ts, Ss);
    X(1,:) = X0;
    R = zeros(Ts, Rs);    
    IntG = zeros(Ts, Rs);
    
    for tIdx = 1:Ts-1
        [X(tIdx+1,:), R(tIdx+1,:), IntG(tIdx+1,:)] = opts.fSimulate(X(tIdx,:), theta, dts(tIdx));
    end
end
