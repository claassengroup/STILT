function [ data ] = stlLoadData( filePath, varargin )
% Loads and validates measurements in Stilt format
% stlLoadData( matFile, matVariableName )
    [~,~,ext] = fileparts(filePath);
    
    switch ext
        case '.mat'
            if (nargin<2)
                error('Name of variable in mat file not specified.');
            end
            fileContent = load(filePath);
            data = fileContent.(varargin{1});
        otherwise
            error('Not supported yet.');
    end
    stlValidateData(data);
end

