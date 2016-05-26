function [ opts ] = stlSetOptSpecies( opts, speciesID, varargin )
% SETSETOPTSPECIES sets the value of a species' field in the options structure.
%
% opts ... stlOptions structure
% speciesID ... ID of the species to change
% varargin ... pairs of fieldName and values to set
%
% See also stlOptions
    assert(~mod(length(varargin),2), 'Odd number of field/value parameters in ''%s''.', speciesID);
    idxSpecies = find(strcmp({opts.Species.ID}, speciesID));
    assert(numel(idxSpecies)==1, 'Unknown (or duplicate) species id ''%s''.', speciesID);
    
    for pIdx = 1:(length(varargin)/2)
        fieldName = varargin{(pIdx-1)*2+1};
        value = varargin{pIdx*2};
        assert(isfield(opts.Species, fieldName), 'Unknown field name ''%s'' for ''%s''.', fieldName, speciesID);
        opts.Species(idxSpecies).(fieldName) = value;
    end
end

