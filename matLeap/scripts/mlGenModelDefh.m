function genModelDefh( filePath, sNr, thNr )
    outStr = sprintf('#ifndef MODEL_DEF_H_\n#define MODEL_DEF_H_\n#define THNR %d\n#define SPNR %d\n#endif\n', thNr, sNr);

    fHandle = fopen(filePath, 'w');
    fwrite(fHandle, outStr);
    fclose(fHandle);
end

