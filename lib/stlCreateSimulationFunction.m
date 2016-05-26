function [ opts ] = stlCreateSimulationFunction( sbmlModel, opts, modelDir, modelName )

    [Stoich, StoichEduct, Props] = stlSBML2StoichProp(sbmlModel);
    
    [Ss, Rs] = size(Stoich);
    
    funcName = [modelName '_getStoichMatrix'];
    fileID = fopen(fullfile(modelDir, [funcName '.m']),'w');
    fprintf(fileID, 'function [S, Sed] = getStoichMatrix()\n\tS=%s;\n\tSed=%s;\nend', mat2str(Stoich), mat2str(StoichEduct));
    fclose(fileID);
    fStoichMatrix = str2func(funcName);
    
    funcName = [modelName '_computePropensities'];
    fileID = fopen(fullfile(modelDir, [funcName '.m']),'w');
    sProps = [num2cell(1:Rs)',Props]';
    sProps = sprintf('\tP(:,%d) = %s;\n', sProps{:});
    fprintf(fileID, 'function [P] = computePropensities(k,X)\n\tP=zeros(size(X,1),%d);\n%s\nend', Rs, sProps);
    fclose(fileID);
    fPropensities = str2func(funcName);

    opts.fSimulate = @(X0, P, Tend) stlTauLeaping(X0, P, Tend, fStoichMatrix, fPropensities, opts.TauLeapingEpsilon, opts.TauLeapingCriticalNumber);
end

