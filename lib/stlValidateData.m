function stlValidateData( data )
% Validates existence of required fields 'time', 'cellNr', 'inspected'
    required = {'time', 'cellNr', 'inspected'};
    exists = ismember(required, fieldnames(data)); 
    if (any(~exists))
        error('Missing field names: %s', strjoin(required(~exists),', '));
    end
end

