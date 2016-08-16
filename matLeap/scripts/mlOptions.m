function [ opts ] = mlOptions( varargin )
        
    % Set default values
    opts = struct();
    opts.DELTA_PE = 0.05;
    opts.N_STIFF = 100;
    opts.NEWTON_EPS = 0.01;
    opts.MAX_ITER = 30;
    opts.SSA_CONDITION_NUMBER = 10;
    opts.CRITICAL_NUMBER = 10;
    opts.USE_IMPLICIT = true;
    opts.EPSILON = 0.03;
    opts.N_C = 10;
    opts.N_SSA_1 = 100;
    opts.N_SSA_2 = 10;
    opts.SSA_ONLY = false;
    opts.SYM_JAC = false;
    opts.USE_RRE = false;
    
    % Override defaults with specified user inputs
    if (~isempty(varargin))
        inputOpts = varargin2struct(varargin{:});
        inputFieldnames = fieldnames(inputOpts);
        isFieldValid = ismember(inputFieldnames, fieldnames(opts));
        if ( any(~isFieldValid) )
            error('Unknown options: %s', strjoin(inputFieldnames{~isFieldValid},','));
        end
        opts=mergestruct(opts, inputOpts);
    end
    
end

