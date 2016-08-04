%% load sbml
modelName = 'BirthDeath';
modelDir = 'models/BirthDeath';

sbmlPath = fullfile(modelDir,[modelName '.xml']);
sbmlModel = TranslateSBML(sbmlPath);
SpeciesNo = length(sbmlModel.species);

%% generate executable
opts = mlOptions();

mexName = ['ml' modelName];
mexFuncPath = fullfile(modelDir, [mexName '.' mexext]);
disp('Building model executable')
mlPrepareModel(modelDir, mexName, sbmlModel, opts);


%% configure run
trajectoryCount = 1;
timeIntervals = 10;

modelStopTime = 3;
X0 = 30;
k_p = 50; % true birth rate
g_p = 5; % true death rate
theta = [k_p g_p];


%% run simulation
disp('Running simulation')
[ X, outTimes, SimTime, r, G, simStats ] = mlSimulate( str2func(mexName), modelStopTime, timeIntervals, trajectoryCount, X0, theta, opts);
disp('Done')

%% show trajectory
figure
plot(outTimes, X)
%% define gamma priors (mean correct 
close all
syms a b


% birth rate
y=solve([a/b==k_p*normrnd(1,0.1), a/b^2==k_p*normrnd(1,0.2)],[a,b]); 
alpha = eval(y.a);
beta = eval(y.b);

% death rate
y2=solve([a/b==g_p*normrnd(1,0.1), a/b^2==g_p*normrnd(1,0.2)],[a,b]);

alpha2 = eval(y2.a);
beta2 = eval(y2.b);


% show posteriors
figure
subplot(121)
hold on
pal = lines(timeIntervals+1);
obsVec = [1,5,11]; % times to show the posterior
for k=obsVec
    i_R = sum(r(1,1:k)); % cumulative r
    i_G = sum(G(1,1:k)); % cumulative G
    [x,y] = fplot(@(x) gampdf(x, alpha+i_R,1/(beta+i_G)), [0,2*k_p], 'LineWidth',2,'Color',pal(k,:));
    plot(x,y, 'LineWidth', 2)
end

h=line(ones(1,2)*k_p, get(gca,'YLim'));
set(h,'Color','k','LineWidth',2)
title('Death rate')
set(gca,'FontSize',14)

subplot(122)
hold on
for k=obsVec
    i_R = sum(r(2,1:k)); % cumulative r
    i_G = sum(G(2,1:k)); % cumulative G
    [x,y]=fplot(@(x) gampdf(x, alpha2+i_R,1/(beta2+i_G)), [0,3*g_p], 'LineWidth',2,'Color',pal(k,:));
    plot(x,y, 'LineWidth', 2)
end
h=line(ones(1,2)*g_p, get(gca,'YLim'));
set(h,'Color','k','LineWidth',2)
ts = linspace(0,modelStopTime,timeIntervals+1);
title('Birth rate')
legend(strsplit(strtrim(sprintf('t=%0.2f\n',ts(obsVec))),'\n'), 'True')
set(gca,'FontSize',14)
