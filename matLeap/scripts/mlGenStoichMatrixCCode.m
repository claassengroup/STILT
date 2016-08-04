function genStoichMatrixCCode(S, Sed, file)

output = sprintf('void computeStoichiometricMatrix(int *S, int *Sed)\n{\n');
ST = S';
for k=1:numel(S)
    output = [output, sprintf('S[%d] = %d;\n', k-1, ST(k))];
end

SedT = Sed';
for k=1:numel(S)
    output = [output, sprintf('Sed[%d] = %d;\n', k-1, SedT(k))];
end

output = [output, '}'];

fh = fopen(file,'w');
fwrite(fh, output);
fclose(fh);
