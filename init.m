
%% Configure setup
% Save the STILT paths which are added by this init script permanently.
SAVE_PATH = false; % do not save STILT to the Matlab working path
% SAVE_PATH = true; % save STILT to the Matlab working path permanently

% Compile Matlab's poissrnd function. Speeds up simulation. Cannot be
% combined with compilation of entire simulation function (stlOptions.Compile)
COMPILE_POISSONRND = true; % compile
% COMPILE_POISSONRND = false; % do not compile

%% Add STILT paths
[InstallDir,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(InstallDir));
if (SAVE_PATH)
    savepath;
end

%% Compile poissrnd
if (COMPILE_POISSONRND)
    PoissrndPath = fullfile(InstallDir,'thirdparty','poissrnd');
    if (~exist(PoissrndPath,'file'))
        cfg = coder.config('mex');
        Lambdatype = coder.newtype('double', [Inf 1], [1 0]);
        eval(['codegen poissrnd -args {Lambdatype} -config cfg -o ' PoissrndPath]);   
    end
end

%% run matLeap init
run(fullfile('matLeap','init.m'))
