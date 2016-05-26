
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

    for k=1:SampleNr
        simtree = stlSimulateTree(X0s(k,:), theta, NGen, obs.time, opts, sbmlModel);
        sim_data{k} = convertSimToData(simtree, obssigmas, opts);
    end

    if mean_type == 1 % geometric
        fMean = @(ll, num) ll / num; 
    else % arithmetic
        fMean = @(ll, num) ll - log(num); 
    end    
    
    %% compute the conditional log likelihoods for each synthetic data set 
    synlogLiks = zeros(SampleNr,1);
    parfor k=1:SampleNr
        tic
        res = stlParticleFilterTree(sbmlModel, sim_data{k}, opts, theta);
        synlogLiks(k) = fMean(res.margLogLik,length(res.meanLogLik)); 
        t=toc;
        fprintf('l=%f [%f s]\n', synlogLiks(k), t)
    end
    
    %% compute the conditional log likelihoods for data
    res = stlParticleFilterTree(sbmlModel, data, opts, theta);
    dataLogLik = fMean(res.margLogLik,length(res.meanLogLik)); 
    
    % Probability data | synTrees
    [f,x]=ecdf(synlogLiks);
    [x, ix] = unique(x);
    f = f(ix);
    pData=interp1(x,f,dataLogLik);
end

function [data] = convertSimToData(simtree, sigmas, opts)

    extractSpecies = find(cellfun(@(c) ~isempty(c), {opts.Species.DataIdx}));
    
    Cs = size(simtree,1);
    data = struct();
    data.time = [];
    data.cellNr = [];
    for sIdx = 1:length(extractSpecies)
        data.(opts.Species(extractSpecies(sIdx)).Name) = [];
        data.([opts.Species(extractSpecies(sIdx)).Name 'Sigma']) = [];
    end    
    for cIdx = 1:Cs
        timeStartIdx = 2-(cIdx==1);
        curTime = simtree{cIdx,2}(timeStartIdx:end);
        data.time = [data.time;curTime]; 
        data.cellNr = [data.cellNr;cIdx*ones(length(curTime),1)]; 
        for sIdx = 1:length(extractSpecies)
        	data.(opts.Species(extractSpecies(sIdx)).Name) = [data.(opts.Species(extractSpecies(sIdx)).Name);simtree{cIdx,3}(timeStartIdx:end,extractSpecies(sIdx))];
            data.([opts.Species(extractSpecies(sIdx)).Name 'Sigma']) = [data.([opts.Species(extractSpecies(sIdx)).Name 'Sigma']);sigmas.(opts.Species(extractSpecies(sIdx)).Name){cIdx}];
        end
    end
    data.inspected = ones(length(data.time),1);
end

