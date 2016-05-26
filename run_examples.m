
% Load SBML model and data
exampleDir = 'examples\Nanog\';
model1Name = 'NoFeedback';
sbmlModel = TranslateSBML(fullfile(exampleDir, [model1Name,'.xml']));
data = stlLoadData(fullfile(exampleDir,[model1Name '_Sim_Protein.mat']),'data');

% Create options structure
opts = stlOptions(sbmlModel, data, 1e4);

% Set parameter priors
opts = stlSetOptParam(opts, 'k1', 'GammaPrior', [5, 2]); % DNA on rate 
opts = stlSetOptParam(opts, 'k2', 'GammaPrior', [5, 2]); % DNA off rate 
opts = stlSetOptParam(opts, 'k3', 'GammaPrior', [12, 0.2]); % mRNA birth rate
opts = stlSetOptParam(opts, 'k4', 'GammaPrior', [7, 1.0]); % mRNA death rate
opts = stlSetOptParam(opts, 'k5', 'GammaPrior', [4.3, 0.005]); % protein birth rate 
opts = stlSetOptParam(opts, 'k6', 'GammaPrior', [5, 5]); % protein death rate

% Configure species' initalization
opts = stlSetOptSpecies(opts, 'DNA_on', 'Init', 'binornd(1, 0.5, [P,1])');
opts = stlSetOptSpecies(opts, 'DNA_off', 'Init', '1-DNA_on');
opts = stlSetOptSpecies(opts, 'RNA', 'Init', 'randi(50,P,1)');
opts = stlSetOptSpecies(opts, 'Protein', 'Init', 'max(0,round(normrnd(Protein1, ProteinSigma1, [P, 1])))');

% Configure cell division
opts = stlSetOptSpecies(opts, 'RNA', 'Division', 'binomial');
opts = stlSetOptSpecies(opts, 'Protein', 'Division', 'binomial');

% Configure plotting for DNA_on/off. Median only. Must be contained in
% opts.LatentQuantiles
opts = stlSetOptSpecies(opts, 'DNA_on', 'PlotQuantiles', 0.5);
opts = stlSetOptSpecies(opts, 'DNA_off', 'PlotQuantiles', 0.5);

% Create forward simulation files
opts = stlCreateSimulationFunction(sbmlModel, opts, exampleDir, model1Name);

% Specify output directory
opts.OutDir = fullfile(exampleDir, [model1Name '_Results']);

% Run particle filter
res1 = stlParticleFilterTree(sbmlModel, data, opts);

% Plot posterior of parameters
stlPlotPosterior(res1.c, sbmlModel, opts, []);

% Plot measurement and prediction for all species
stlPlotTimecourseSpeciesQuantilesTree(sbmlModel, res1.obs, opts, res1.Q, []);

%% Negative Feedback (on no-feedback data)

% Load SBML model and data
exampleDir = 'examples\Nanog\';
model2Name = 'NegFeedback';
sbmlModel = TranslateSBML(fullfile(exampleDir, [model2Name,'.xml']));
data = stlLoadData(fullfile(exampleDir,[model1Name '_Sim_Protein.mat']),'data');

% Create options structure
opts = stlOptions(sbmlModel, data, 1e4);

% Set parameter priors
opts = stlSetOptParam(opts, 'k1', 'GammaPrior', [5, 0.8]); % DNA on rate 
opts = stlSetOptParam(opts, 'k2', 'GammaPrior', [7, 1e5]); % DNA off rate 
opts = stlSetOptParam(opts, 'k3', 'GammaPrior', [12, 0.2]); % mRNA birth rate
opts = stlSetOptParam(opts, 'k4', 'GammaPrior', [7, 1.0]); % mRNA death rate
opts = stlSetOptParam(opts, 'k5', 'GammaPrior', [4.3, 0.005]); % protein birth rate 
opts = stlSetOptParam(opts, 'k6', 'GammaPrior', [5, 5]); % protein death rate

% Configure species' initalization
opts = stlSetOptSpecies(opts, 'DNA_on', 'Init', 'binornd(1, 0.5, [P,1])');
opts = stlSetOptSpecies(opts, 'DNA_off', 'Init', '1-DNA_on');
opts = stlSetOptSpecies(opts, 'RNA', 'Init', 'randi(50,P,1)');
opts = stlSetOptSpecies(opts, 'Protein', 'Init', 'max(0,round(normrnd(Protein1, ProteinSigma1, [P, 1])))');

% Configure cell division
opts = stlSetOptSpecies(opts, 'RNA', 'Division', 'binomial');
opts = stlSetOptSpecies(opts, 'Protein', 'Division', 'binomial');

% Configure plotting for DNA_on/off. Median only. Must be contained in
% opts.LatentQuantiles
opts = stlSetOptSpecies(opts, 'DNA_on', 'PlotQuantiles', 0.5);
opts = stlSetOptSpecies(opts, 'DNA_off', 'PlotQuantiles', 0.5);

% Create forward simulation files
opts = stlCreateSimulationFunction(sbmlModel, opts, exampleDir, model1Name);

% Specify output directory
opts.OutDir = fullfile(exampleDir, [model2Name '_Results']);

% Run particle filter
res2 = stlParticleFilterTree(sbmlModel, data, opts);

% Plot posterior of parameters
stlPlotPosterior(res2.c, sbmlModel, opts, []);

% Plot measurement and prediction for all species
stlPlotTimecourseSpeciesQuantilesTree(sbmlModel, res2.obs, opts, res2.Q, []);

%% Model comparison
bf = res1.margLogLik - res2.margLogLik;
fprintf('log(Bayes factor) [True over False] = %d\n', bf);


