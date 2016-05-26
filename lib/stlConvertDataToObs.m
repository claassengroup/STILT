function [ obs, sigmas ] = stlConvertDataToObs( data, opts )
% STLCONVERTDATATOOBS structures the list of observations (data) by cell

    uniqueTPs = unique(data.time);
    CNr = length(unique(data.cellNr));
    cellNrs = unique([data(:).cellNr]);
    dataSpecies = find(cellfun(@(c) ~isempty(c), {opts.Species.DataIdx}));
    
    obs = struct();
    sigmas = struct();
    for k=1:CNr
        ix = cellNrs(k);
        tmp = find(data.cellNr==ix);
        obs.cellNr(k) = ix;
        [~,tpIdcs]=ismember(data.time(tmp), uniqueTPs);
        if (any(diff(tpIdcs)~=1))
            error('Missing timepoint in cell %d.', k);
        end
        obs.time{k} = reshape(data.time(tmp),[],1);
        for s = 1:length(dataSpecies)
            speciesName=opts.Species(dataSpecies(s)).Name;
            obs.(speciesName){k} =  data.(speciesName)(tmp);
            speciesSigma = data.([speciesName 'Sigma']);
            if (numel(speciesSigma) == numel(data.(speciesName)))
                obs.([speciesName 'Sigma']){k} =  speciesSigma(tmp);
            else
                obs.([speciesName 'Sigma']){k} = speciesSigma * ones(size(tmp));
            end
            sigmas.(speciesName){k} = obs.([speciesName 'Sigma']){k};
        end
        if ~isfield(data, 'inspected')
            obs.inspected{k} = ones(size(tmp));
        else
            obs.inspected{k} = data.inspected(tmp);
        end
    end
end

