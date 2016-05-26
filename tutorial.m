%% STILT 1.0: Installation and tutorial

%% Introduction
% This document contains information about the installation and usage of
% STILT.  For details on the method please see
% TODO: citation

init
%% Installation
% The STILT toolbox is a set of Matlab scripts/functions and is 'installed'
% by copying and adding the directory to Matlab's search path. Furthermore,
% the stochastic simulation within STILT can (optionally) be accelerated by compiling 
% Matlab's Poisson random number function (poissrnd). Installation and
% compilation can be performed by running the script _init.m_ in STILTs top-level directory.
% 
% Before you run init.m, you can configure the following variables (found in init.m):
%
% * SAVE_PATH = true | false 
% * COMPILE_POISSONRND = true | false  
%
% STILT uses libSBML-5.12.0 to load SBML models. If another compatible version of this
% library is installed on the Matlab path, the directory /thirdparty/libSBML-5.12.0-matlab can
% be deleted.
%

%% SBML model
%
% STILT reads SBML files via
% <http://sbml.org/Software/libSBML/docs/matlab-api/ libSBML> and uses the
% species, parameters and reactions definitions contained within. The information is used for 
% STILT's options initialization and for the generation of a model specific 
% forward simulation function (from stoichiometry and reaction definition). 
%
% Biochemical species for which there is no net change in copy number should
% nevertheless be specified as educt and product when necessary for a reaction to proceed. 
% For example, for the reaction *X + Y --> X + Z*
% the species *X* is both educt and product.
% The definition of reaction rates is imported from _kineticlaw.formula_ field
% of the reaction definition with in the SBML file. Use Latex syntax in the names of
% species and parameters for plotting.  See examples/Nanog/NoFeedback.xml
% and examples/Nanog/NegativeFeedback.xml for example SBML model configuration.

%% Data
% 
% Experimental data used for parameter inference is specified as a Matlab structure with 
% the fields _time_, _cellNr_ and _inspected_, which must be vectors of equal length.
% The _ith_ entry of each of these vectors specifies a single observation at time _time(i)_
% of the cell with index _cellNr(i)_.  
% Cellular genealogy is introduced via the field _cellNr_, such that 
% the parent cell with index _j_ gives rise to daughter cells with indices _2*j_ and
% _2*j+1_, e.g. cell 1 has daughter cells 2 and 3. 
% If the observation at time _i_ is to be censored (i.e. not used for inference),
% the field _inspected(i)_ should be set to 0, otherwise 1. 
% For each observed species one field with the species' name and one with the
% species' and the suffix _Sigma_ containing the measurement noise magnitudes at each time point
% should also be included in the data structure (see example below). 
% Example structure with a single observed species named "Protein" and 147 observations:
%
%            time: [1x147 double]
%          cellNr: [1x147 double]
%         Protein: [1x147 double]
%       inspected: [1x147 double]
%    ProteinSigma: 200

%% Configuration
% The tutorial contains the steps necessary to perform model selection
% on linage tree data as shown in the publication. 

%%
% <html>
% <h3>Working directory, model, data</h3>
% </html>
%
% STILT requires an SBML model definition and the measurement data as input.
% Data an either be specified by creating a data structure with the fields above,
% or by using the function _stlLoadData_ to load a .mat file containing a structure named _data_ with the 
% required format. _stlLoadData_ also verifies that the data are correctly
% formatted.

% Our working directory. Everything (including results) stored in there.
exampleDir = 'examples\Nanog\';
% Load SBML via libSBML
sbmlModelNoFeedback = TranslateSBML(fullfile(exampleDir, 'NoFeedback.xml'));
% Load data from mat file.
% data = stlLoadData(fullfile(exampleDir,[modelName '_Sim_Protein.mat']),'data');

%%
% <html>
% <h3>Create options structure</h3>
% </html>
%
% The function _stlOptions_ creates a structure with all options available
% in STILT set to their default values.  It requires the following three arguments
% in the specified order: i) the _sbmlModel_, ii) the _data_ structure, and iii) the number of 
% particles to use for inference.

% Create options template from model, data and particle number
optsNoFeedback = stlOptions(sbmlModelNoFeedback, [], 7e5);
% Set output directory for results (figures, data, etc will be in subdirectories)
optsNoFeedback.OutDir = fullfile(exampleDir, 'NoFeedbackNull_Results');
%%
% <html>
% <h3>Configure species' initalization</h3>
% </html>
%
% The function _stlSetOptSpecies_(_opt_, _ParameterName_, _FieldName_, _FieldValue_) can be
% used to configure the behavior of individual species in the options
% structure. Options are specified as name/value (_FieldName_,_FieldValue_) pairs.
% The field 'Init' contains a Matlab command (as string), which
% is evaluated at runtime for particle initialzation. The number of particles
% specified at run time can be accessed with the variable _P_. 
% The variables '<SpeciesName>1' and
% '<SpeciesName>Sigma1' contain the first data point for each measured
% species. Previously initialized variables (each a P-by-1 vector) can be
% referenced via their species' names.

optsNoFeedback = stlSetOptSpecies(optsNoFeedback, 'DNA_on', 'Init', 'binornd(1, 0.5, [P,1])');
optsNoFeedback = stlSetOptSpecies(optsNoFeedback, 'DNA_off', 'Init', '1-DNA_on');
optsNoFeedback = stlSetOptSpecies(optsNoFeedback, 'RNA', 'Init', 'randi(50,P,1)');
optsNoFeedback = stlSetOptSpecies(optsNoFeedback, 'Protein', 'Init', 'max(0,round(normrnd(Protein1, ProteinSigma1, [P, 1])))');

%%
% <html>
% <h3>Configure cell division</h3>
% </html>
%
% STILT currently supports two modes of cell division possible for each simulated species.
% They may be configured for each species individually via the field 'Division' in the STILT options structure: 
%
% * 'copy' ... (default) the mother's content is copied to both daughter cells (e.g. for DNA)
% * 'binomial' ... the mother's content M is split binomially between the 
% daughters. B(M, 0.5) (if M > 100 approximated with a normal distribution.)
% * 'binomial' with the oarameter _p_ specified using the additional
% name/value pair _DivisionParameters_.
optsNoFeedback = stlSetOptSpecies(optsNoFeedback, 'RNA', 'Division', 'binomial');
% opts_NoFeedback = stlSetOptSpecies(opts_NoFeedback, 'Protein', 'Division', 'binomial'); % default p = 0.5
optsNoFeedback = stlSetOptSpecies(optsNoFeedback, 'Protein', 'Division', 'binomial', 'DivisionParameters', 'p=0.5'); % specify p explicitly

%%
% <html>
% <h3>Prepare forward simulation</h3>
% </html>
%
% STILT's particle filter is based on stochastic simulation (performed by
% the function _stlTauLeaping_). The function _stlCreateSimulationFunction_ generates code for performing 
% stochastic simulation using the information contained in _sbmlModel_ (i.e. the stoichiometry and reaction propensity).
% The field_opts.fSimulate_ contains a function handle to these simulation files.

optsNoFeedback = stlCreateSimulationFunction(sbmlModelNoFeedback, optsNoFeedback, exampleDir, 'NoFeedback');

%% Generate synthetic data using the No Feedback model
% initial conditions (DNA on, DNA off, mRNA, protein)
X0 = [0, 1, 50, 1e4];
% model parameters (k_on, k_off, k_m, g_m, k_p, g_p)
theta = [1, 1, 40, 3, 300, 0.4];
% observe only the protein
optsNoFeedback = stlSetOptSpecies(optsNoFeedback, 'Protein', 'Observed', true);
optsNoFeedback = stlSetOptSpecies(optsNoFeedback, 'Protein', 'DataIdx', 1);

timeParams = {12,1,0.5}; % mean, std of cell life times and the dt parameter for observation frequency
simData_NoFeedback = stlSimulateDataSet(X0, theta, 3, {12,1,0.5}, optsNoFeedback, sbmlModelNoFeedback, 500);

%%
% <html>
% <h3>Configure parameter priors</h3>
% </html>
%
% The function _stlOptions_ creates empty entries for each of the parameters in
% _sbmlModel_. Set individual parameter's fields using the function 
% _stlSetOptParam_(_opt_, _ParameterName_, _FieldName_, _FieldValue_) as shown below. 
% The hyperparameters of the gamma prior distributions ($$ \alpha $$, $$ \beta $$)
% for each model parameter should be configured as shown below.

optsNoFeedback = stlSetOptParam(optsNoFeedback, 'k_on', 'GammaPrior', [5, 2]); % DNA on rate 
optsNoFeedback = stlSetOptParam(optsNoFeedback, 'k_off', 'GammaPrior', [5, 2]); % DNA off rate 
optsNoFeedback = stlSetOptParam(optsNoFeedback, 'k_m', 'GammaPrior', [12, 0.2]); % mRNA birth rate
optsNoFeedback = stlSetOptParam(optsNoFeedback, 'g_m', 'GammaPrior', [7, 1.0]); % mRNA death rate
optsNoFeedback = stlSetOptParam(optsNoFeedback, 'k_p', 'GammaPrior', [4.3, 0.005]); % protein birth rate 
optsNoFeedback = stlSetOptParam(optsNoFeedback, 'g_p', 'GammaPrior', [5, 5]); % protein death rate

% Set the field _TrueParamValuesInModel_ in the options for STILT to _true_
% if the numeric parameter values from SBML are assumed to be correct
% (i.e. when using synthetic data for which the correct parameters are known)
% and should be plotted, e.g. in the parameter posterior plot.  




%% Change the parameters _ShowPlotsWhileFiltering_ and _SavePlotsWhileFiltering_ to
% generate and save plots, respectively, while running the algorithm.
optsNoFeedback.ShowPlotsWhileFiltering = true;
optsNoFeedback.SavePlotsWhileFiltering = true;

%% Running the particle filter
% The function _stlParticleFilterTree_ performs the parameter inference.
% As a result, it generates a results structure _res_ which consists of the
% following fields:
%
% * _c_ ... a matrix of samples from the parameters' posterior
% distributions (particles-by-parameters)
% * _w_ ... a vector of weights for the individual particles
% (particles-by-1)
% * _alpha_, _beta_  ... the hyperparameters of each particle's gamma distribution, (particles-by-1) each
% * _Q_ ... a matrix of quantiles of the distributions of the simulated
% system state for each cell (timepoints-by-species-by-quantiles)
% * _obs_ ... _data_ re-structured cell oriented for plotting


% configure _LatentQuantilesPerTimestep_ >= 1 to configure the number of
% timesteps to use for interpolating between observations (helps smooth
% estimated trajectories).
opts.LatentQuantilesPerTimestep=1;
% run the particle filter with the No Feedback model
optsNoFeedback.ParticleNr = 7e5;
resNoFeedback_Null = stlParticleFilterTree(sbmlModelNoFeedback, simData_NoFeedback, optsNoFeedback);

%% Plot

%%
% <html>
% <h3>Posterior of parameters</h3>
% </html>
%
% _stlPlotPosterior_ plots samples from each parameter's posterior (black),
% its true value (if known, red vertical line), the prior (red dashed line)
% and the prior mean (blue vertical line).

% resample to smooth the histograms
stlPlotPosterior(resNoFeedback_Null.c, sbmlModelNoFeedback, optsNoFeedback, []);

%%
% 
% <<parameter_posterior.png>>
%

%%
% <html>
% <h3>Plot prediction and (if available) measurements for all species</h3>
% </html>
%
% _stlPlotTimecourseSpeciesQuantilesTree_ generates a plot with the time course 
% of simulation predictions (and measurements) arranged per cell in a tree
% format for each individual species.

% stlPlotTimecourseSpeciesQuantilesTree(sbmlModelNoFeedback, resNoFeedback_Null.obs, opts_NoFeedback, resNoFeedback_Null.Q, []);

%%
% 
% <<speciesquantiletree.png>>
%

%% Fit the data from the no feedback model with the negative feedback model
% Our working directory. Everything (including results) stored in there.
% Load SBML via libSBML
sbmlModelNegativeFeedback = TranslateSBML(fullfile(exampleDir, 'NegFeedback.xml'));

% Create options template from model, data and particle number
optsNegFeedback = stlOptions(sbmlModelNegativeFeedback, simData_NoFeedback, 7e5);
% Set output directory for results (figures, data, etc will be in subdirectories)
optsNegFeedback.OutDir = fullfile(exampleDir, 'NoFeedbackNeg_Results');
%%
% set parameter prior parameters
optsNegFeedback = stlSetOptParam(optsNegFeedback, 'k_on', 'GammaPrior', [5, 0.8]); % DNA on rate 
optsNegFeedback = stlSetOptParam(optsNegFeedback, 'k_off', 'GammaPrior', [7, 1e5]); % DNA off rate 
optsNegFeedback = stlSetOptParam(optsNegFeedback, 'k_m', 'GammaPrior', [12, 0.2]); % mRNA birth rate
optsNegFeedback = stlSetOptParam(optsNegFeedback, 'g_m', 'GammaPrior', [7, 1.0]); % mRNA death rate
optsNegFeedback = stlSetOptParam(optsNegFeedback, 'k_p', 'GammaPrior', [4.3, 0.005]); % protein birth rate 
optsNegFeedback = stlSetOptParam(optsNegFeedback, 'g_p', 'GammaPrior', [5, 5]); % protein death rate

% set cell division parameters
optsNegFeedback = stlSetOptSpecies(optsNegFeedback, 'DNA_on', 'Init', 'binornd(1, 0.5, [P,1])');
optsNegFeedback = stlSetOptSpecies(optsNegFeedback, 'DNA_off', 'Init', '1-DNA_on');
optsNegFeedback = stlSetOptSpecies(optsNegFeedback, 'RNA', 'Init', 'randi(50,P,1)');
optsNegFeedback = stlSetOptSpecies(optsNegFeedback, 'Protein', 'Init', 'max(0,round(normrnd(Protein1, ProteinSigma1, [P, 1])))');

% create simulation functions
optsNegFeedback = stlCreateSimulationFunction(sbmlModelNegativeFeedback, optsNegFeedback, exampleDir, 'NegFeedback');

% Change the parameters _ShowPlotsWhileFiltering_ and _SavePlotsWhileFiltering_ to
% generate and save plots, respectively, while running the algorithm.
optsNegFeedback.ShowPlotsWhileFiltering = 1;
optsNegFeedback.SavePlotsWhileFiltering = 1;

%% run the particle filter with the Negative Feedback model and the No Feedback data
optsNegFeedback.ParticleNr = 7e5;
resNoFeedback_Neg = stlParticleFilterTree(sbmlModelNegativeFeedback, simData_NoFeedback, optsNegFeedback);

%% Generate Data from the Negative Feedback model

% initial conditions
X0 = [0, 1, 50, 3e3];
% model parameters (k_on, k_off, k_m, g_m, k_p, g_p)
theta = [10, 1e-4, 100, 3, 1000, 0.4];

timeParams = {12,1,0.5}; % mean, std of cell life times and the dt parameter for observation frequency
simData_NegativeFeedback = stlSimulateDataSet(X0, theta, 3, timeParams, optsNegFeedback, sbmlModelNegativeFeedback, 500);

%% run the particle filter with the Negative Feedback model and Neg Feedback data
optsNegFeedback.ParticleNr = 7e5;
optsNegFeedback.OutDir = fullfile(exampleDir, 'NegFeedbackNeg_Results');
resNegFeedback_Neg = stlParticleFilterTree(sbmlModelNegativeFeedback, simData_NegativeFeedback, optsNegFeedback);

%% fit the Negative Feedback data with the No Feedback model

% Set output directory for results (figures, data, etc will be in subdirectories)
optsNoFeedback.OutDir = fullfile(exampleDir, 'NegFeedbackNull_Results');
resNegFeedback_Null = stlParticleFilterTree(sbmlModelNoFeedback, simData_NegativeFeedback, optsNoFeedback);

%% Model comparison with Bayes factors
% 
% Perform Bayesian model comparison by computing the Bayes factor between two
% models directly from the marginal log likelihood of the data given the
% model P(D_i|M_i) returned in _stlParticleFilterTree's_ results.

disp('No feedback model is true:')
fprintf('Log P(no feedback) - Log P(neg feedback): %f\n', resNoFeedback_Null.margLogLik - resNoFeedback_Neg.margLogLik)

disp('Negative feedback model is true:')
fprintf('Log P(neg feedback) - Log P(no feedback): %f\n', resNegFeedback_Neg.margLogLik - resNegFeedback_Null.margLogLik)
