function data = simulateDataSet(X0, theta, NGen, timesParams, opts, sbmlModel, sigmas)

% generate an artificial data set using the specified model
%
% X0 ... Initial state of initial cell
% theta ... vector of parameters
% NGen ... number of generations
% timeParams ... cell array (1x3) giving mean and std of life times of each
%  cell, and the observation interval dt.  If timeParams is a vector then
%  use this vector of times instead.
% opts ... STILT options structure
% sbmlModel ... SBML model
% sigmas... a vector of constant measurement errors to add to each
% observation
% Returns:
% data set containing the fields:
%  - time
%  - cellNr
%  - <SpeciesName> for each observable species
%  - <SpeciesName>Sigma for measurement error of each obs species
%  - inspected: a flag indicating whether the data are observed

%% generate a cell array of times for each of the cells

Cs = 2^NGen-1;
mu_death = timesParams{1};
sig_death = timesParams{2};
dt = timesParams{3};

deathTime{1} = max(1, normrnd(mu_death, sig_death));
times{1} = reshape(0:dt:(dt*floor(deathTime{1}/dt)),[],1);
for ixCell = 2:Cs
    motherCell = floor(ixCell/2);
    deathTime{ixCell} = deathTime{motherCell} + max(1, normrnd(mu_death, sig_death));
    times{ixCell} = reshape((times{motherCell}(end) + dt):dt:(deathTime{ixCell}),[],1);
end
%%
% simulate the tree
[ Cells ] = stlSimulateTree(X0, theta, NGen, times, opts, sbmlModel);

observedSpecies = find([opts.Species.Observed]);
% create fields

nObs = length(observedSpecies);
for i=1:nObs
    ixObs = observedSpecies(i);
    data.(opts.Species(ixObs).ID) = [];
    data.(strcat(opts.Species(ixObs).ID, 'Sigma')) = [];
end
data.time = [];
data.cellNr = [];

% loop over cells

for ixCell = 1:Cs
    X = Cells{ixCell,3}; % states
    nRows = size(X,1);
    for i=1:nObs
        ixObs = observedSpecies(i);
        data.(opts.Species(ixObs).ID) = [data.(opts.Species(ixObs).ID); X(:,ixObs) + ...
            normrnd(0, sigmas(i), size(X(:,1))) ];
        data.(strcat(opts.Species(ixObs).ID, 'Sigma')) = [data.(strcat(opts.Species(ixObs).ID, 'Sigma')); sigmas(i)*ones(nRows,1)];
    end
    data.time = [data.time; Cells{ixCell, 2}];
    data.cellNr = [data.cellNr; ixCell*ones(nRows,1)];
end

data.inspected = ones(size(data.time)); % all inspected

        
    
