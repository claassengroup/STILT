function [ opts ] = stlSetOptParam( opts, parameterID, fieldName, value )
% SETSETOPTPARAMETER sets the value of a species' field in the options structure.
%
% opts ... stlOptions structure
% parameterID ... ID of the parameter to change, string
% fieldName ... Name of the option field within parameter
% value ... Value to be set
%
% See also stlOptions

    idxParam = find(strcmp({opts.Parameters.ID}, parameterID));
    assert(numel(idxParam)==1, 'Unknown (or duplicate) parameter id ''%s''.', parameterID);
    assert(isfield(opts.Parameters, fieldName), 'Unknown field name ''%s'' for ''%s''.', fieldName, parameterID);
    opts.Parameters(idxParam).(fieldName) = value;
end

