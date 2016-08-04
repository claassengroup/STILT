
%% Configure setup
% Save the matLeap paths which are added by this init script permanently.
SAVE_PATH = 0;

%% Add STILT paths
[InstallDir,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(InstallDir));
if (SAVE_PATH)
    savepath;
end
