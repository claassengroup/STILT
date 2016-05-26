function data = simulateDataSet(X0, theta, NGen, times, opts, sbmlModel, sigmas)

% generate an artificial data set using the specified model
%
% X0 ... Initial state of initial cell
% theta ... vector of parameters
% NGen ... number of generations
% times ... cell array (cells-by-1) of time-vectors per cell
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
Cs = 2^NGen-1;

for ixCell = 1:Cs
    X = Cells{ixCell,3}; % states
    nRows = size(X,1);
    times = Cells{ixCell,2};
    for i=1:nObs
        ixObs = observedSpecies(i);
        data.(opts.Species(ixObs).ID) = [data.(opts.Species(ixObs).ID); X(:,ixObs)];
        data.(opts.Species(ixObs).ID) = [data.(opts.Species(ixObs).ID); sigmas(i)*ones(nRows,1)];
    end
    data.time = [data.time; Cells{ixCell, 2}];
    data.cellNr = [data.cellNr; ixCell*ones(nRows,1)];
end

data.inspected = ones(size(data.time)); % all inspected

        
    
