%% Configure setup
% Save changes to path and startup.m
SAVE_PERMANENTLY = 1 ; % 1 to save
%% Add paths
[InstallDir,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(InstallDir));
if (SAVE_PERMANENTLY)
    savepath;
end
%% Set environment variable with install directory
setenv('matLeap1_InstallDir',InstallDir);
if (SAVE_PERMANENTLY)
    sDir = fileparts(which('startup'));
    sPath = fullfile(sDir, 'startup.m');
    sContent = fileread(sPath);
    if (isempty(regexp(sContent, 'matLeap1_InstallDir')))
        fileID = fopen(sPath,'a');
        fprintf(fileID,'\nsetenv(''matLeap1_InstallDir'',''%s'');\n',InstallDir);
        fclose(fileID);
    end
end