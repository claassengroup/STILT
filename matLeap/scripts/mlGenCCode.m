function genCCode(name, snippet, inputVars, outputVar, dims, constIn, outDoublePointer, isPropFunc, file)

    inStr = '';
    assignInStr = '';

    constInStr = '';
    if (constIn)
        constInStr = 'const ';
    end
    
    % one-off for the scalar case
    if prod(dims)==1        
        snippet = regexprep(snippet, '(\w+) = (.*)', [outputVar '[0][0] = $2']);
    end

    for k=1:length(inputVars)
        if k ~= length(inputVars)
            endTok = ',';
        else
            endTok = '';
        end    
        inStr = [inStr, sprintf('%sdouble *in%d%s ', constInStr, k, endTok)];
        for l = 1:length(inputVars{k})
            assignInStr = [assignInStr, sprintf('double %s = in%d[%d];\n', inputVars{k}{l}, k, l-1)];
        end
    end
    if (isPropFunc)
        inStr = [inStr, sprintf(',%s int lastReacIdx ', constInStr)];
    end

    if (outDoublePointer)
        outStr = sprintf('double **p_%s', outputVar);
    else
        outStr = sprintf('double *p_%s', outputVar);
    end

    header = sprintf('void %s(%s, %s) {\n', name, inStr, outStr);
    output = sprintf('#include <cstring>\n');
    output = [output, header, assignInStr];

    if (~outDoublePointer)
        if length(dims) == 1
            output = [output, sprintf('double %s[%d];\n', outputVar, dims)];
        else
            output = [output, sprintf('double %s[%d][%d];\n', outputVar, dims(:))];
        end
        output = [output, sprintf('memset(%s, 0.0, sizeof(double)*%d);\n', outputVar, prod(dims))];
    end
    output = [output, snippet, sprintf('\n')];

    if (~outDoublePointer)
        output = [output, sprintf('memcpy(p_%s, &%s[0][0], sizeof(double)*%d);\n', outputVar, outputVar, prod(dims))];
    end
    output = [output, sprintf('\n}')];


    fHandle = fopen(file, 'w');
    fwrite(fHandle, output);
    fclose(fHandle);

end