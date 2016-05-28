
function [pData, synlogLiks, dataLogLik] = stlGoodnessOfFitTest(sbmlModel, data, res, opts)

    mean_type = 2;
    
    %% validate res
    % res can be either the results of the particle filter
    % or a user-specified value for each of the model parameters
    [ffound] = ismember({sbmlModel.parameter.id}, fields(res));
    isResultsStruct = any(~ffound);

    if isResultsStruct
        if ~(isfield(res, 'c') && isfield(res, 'w'))
            error('res must be either result of stlParticleFilter or structure specifying each parameter value')
        end
    end

    %% check for parameters related to goodness-of-fit test
    if isfield(opts, 'GOF_ParticleNr')
        opts.ParticleNr = opts.GOF_ParticleNr;
    else
        opts.ParticleNr = 500;
    end

    if isfield(opts, 'GOF_SampleNr')
        SampleNr = opts.GOF_SampleNr;
    else
        SampleNr = 50;
    end

    %% set theta
    if isResultsStruct
        samples = discretesample(res.w/sum(res.w), length(res.w));
        if ~isfield(opts, 'GOF_PointEstimate') || strcmpi( opts.GOF_PointEstimate, 'median')
            theta = median(res.c(samples,:));
        elseif strcmpi( opts.GOF_PointEstimate, 'mean')
            theta = mean(res.c(samples,:));
        else
            error('Invalid point estimate method. Must be either "median" or "mean".')
        end
    else
        theta = cellfun(@(id) res.(id), {opts.Parameters.ID});
    end

    %% get number of generations from the data
    cells = unique(data.cellNr);
    maxGen = max(floor(log2(cells)));
    minGen = min(floor(log2(cells)));
    NGen = 1 + maxGen - minGen;

    %% convert data for sampling the initial conditions
    observedSpecies = find(cellfun(@(c) c==true, {opts.Species.Observed}));
    [obs, obssigmas] = stlConvertDataToObs(data, opts);

    % extract first time point information for first cell
    cell_ix = find(obs.cellNr == cells(1));
    
    % observed intensities for this cell
    observed1 = zeros(1,length(observedSpecies));
    sigma1 = zeros(1,length(observedSpecies));
    for sIdx=1:length(observedSpecies)
        observed1(sIdx) = obs.(opts.Species(observedSpecies(sIdx)).Name){cell_ix}(1);
        sigma1(sIdx) = obs.([opts.Species(observedSpecies(sIdx)).Name 'Sigma']){cell_ix}(1);
    end
    X0s = stlInitParticles(SampleNr, observed1, sigma1, opts);
    
    %% generate synthetic trees
    sim_data = cell(SampleNr,1);
    disp('Generating synthetic trees...')
    parfor k=1:SampleNr
%         sim_data{k} = stlSimulateDataSet(X0s(k,:), theta, NGen, obs.time, opts, sbmlModel);
        simtree = stlSimulateTree(X0s(k,:), theta, NGen, obs.time, opts, sbmlModel);
        sim_data{k} = stlConvertSimToData(simtree, obssigmas, opts);
    end

    if mean_type == 1 % geometric
        fMean = @(ll, num) ll / num; 
    else % arithmetic
        fMean = @(ll, num) ll - log(num); 
    end    
    
    %% compute the conditional log likelihoods for each synthetic data set 
    synlogLiks = zeros(SampleNr,1);
    
    % disable parallel workers in the tree likelihood computation
    opts.ParallelWorkers = 1;
    disp('Computing likelihood of simulated data sets...')
    tic
    parfor k=1:SampleNr
%         tic
%         fprintf('.')
        res = stlParticleFilterTree(sbmlModel, sim_data{k}, opts, theta);
        synlogLiks(k) = fMean(res.margLogLik,length(res.meanLogLik)); 
%         t=toc;
%         fprintf('x')
%         fprintf('l=%f [%f s]\n', synlogLiks(k), t)
    end
    toc
    %% compute the conditional log likelihoods for data
    res = stlParticleFilterTree(sbmlModel, data, opts, theta);
    dataLogLik = fMean(res.margLogLik,length(res.meanLogLik)); 
    
    % Probability data | synTrees
    [f,x]=ecdf(synlogLiks);
    [x, ix] = unique(x);
    f = f(ix);
    pData=interp1(x,f,dataLogLik);
    if all(dataLogLik<x)
        pData=0;
    elseif all(dataLogLik>x)
        pData=1;
    end
end


