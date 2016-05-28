function [data] = convertSimToData(simtree, sigmas, opts, varargin)

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
        	data.(opts.Species(extractSpecies(sIdx)).Name) = [data.(opts.Species(extractSpecies(sIdx)).Name); ...
                simtree{cIdx,3}(timeStartIdx:end,extractSpecies(sIdx)) + normrnd(0, sigmas.(opts.Species(extractSpecies(sIdx)).Name){cIdx}, size(simtree{cIdx,3}(timeStartIdx:end,extractSpecies(sIdx))))  ];
            data.([opts.Species(extractSpecies(sIdx)).Name 'Sigma']) = [data.([opts.Species(extractSpecies(sIdx)).Name 'Sigma']);sigmas.(opts.Species(extractSpecies(sIdx)).Name){cIdx}];
        end
    end
    data.inspected = ones(length(data.time),1);
end
