function [ opt ] = stlOptions( model, data, particlenr )
% STLOPTIONS returns an options structure template for a specific model and
% data set 
%
% See also stlSetOptParam, stlSetOptSpecies

    opt = struct();

    % Number of Particles
    opt.ParticleNr = particlenr;

    % Copy parameters from sbmlModel
    % GammaPrior ... [alpha, 1/beta] of parameter's gamma distributed prior
    opt.ParamNr = length(model.parameter);
    opt.Parameters = struct('ID', {}, 'Name', {}, 'GammaPrior', {});
    for pIdx = 1:opt.ParamNr
        opt.Parameters(pIdx).ID = model.parameter(pIdx).id;
        opt.Parameters(pIdx).Name = model.parameter(pIdx).name;
    end
    
    % Copy parameters from sbmlModel and set defaults
    % Division ... copy|binomial
    % Init ... Matlab command run at initialization of particles
    % DataIdx ... Number of observed species (set automatically)
    if isempty(data)
        warning('Data missing. ''Species'' options not pre-defined.');
    else
        dataFieldnames = fieldnames(data);
    end
    opt.Species = struct('ID', {}, 'Name', {}, 'Division', {}, 'DivisionParameters', {}, 'Init', {}, 'PlotQuantiles', {}, 'Observed', {}, 'DataIdx', {});
    
    for sIdx = 1:length(model.species)
        opt.Species(sIdx).ID = model.species(sIdx).id;
        opt.Species(sIdx).Name = model.species(sIdx).name;
        opt.Species(sIdx).Division = 'copy';
        opt.Species(sIdx).DivisionParameters = '';
        if ~isempty(data)
            [obs,idx] = ismember(model.species(sIdx).name, dataFieldnames);
            if (obs)
                opt.Species(sIdx).DataIdx = idx;
            end
            opt.Species(sIdx).Observed = obs;
        else
            opt.Species(sIdx).Observed = false;
        end
        
    end

    
    % Output directory
    opt.OutDir = '.';
    
    % Show plots during inference
    opt.ShowPlots = true;
    opt.ShowPlotsWhileFiltering = false;
    
    % Save plots during inference
    opt.SavePlots = true;
    opt.SavePlotsWhileFiltering = false;
    
    % Number of rows / columns for posteror plot
	opt.SubplotRows = ceil(sqrt(length(model.parameter)));
    opt.SubplotColumns = ceil(length(model.parameter)/ceil(sqrt(length(model.parameter))));

    % Save latent states' quantiles on a grid of T time points per cell
    opt.LatentQuantiles = [0.025 0.25 0.50 0.75 0.975];
    opt.LatentQuantilesPerTimestep = 1;
    
    % Indicate whether the SBML files contains the true parameter values (for plotting)
    opt.TrueParamValuesInModel = true;
    
    % Save particles' parameter, weights
    opt.Save = true;
    
    % Save statistics of the accepted particles latent trajectories
    opt.SaveX = false;    
    
    % Function for initialisation of particles
    opt.fInitParticleSpecies = @stlInitParticles;
    
    % TauLeaping parameters
    opt.TauLeapingEpsilon = 0.03;
    opt.TauLeapingCriticalNumber = 10;
    
    % Execute in parallel with Matlab's parfor
    opt.Parallel = false;
    
    % Number of parallel workers (determined automatically if [])
    opt.ParallelWorkers = [];
    
end

