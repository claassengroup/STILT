function [ Trajectories, timeOut, R, IntG, simStats ] = mlSimulate( fTauleaping, modelStopTime, timeIntervals, x0, theta, opts, varargin)
    
if ~isempty(varargin)
    % multiple trajectories
    trajectoryCount = varargin{1};
    
    if ~isrow(theta) || ~isrow(x0)
        error('Cannot specify both matrix X0 / theta AND multiple realizations of the trajectory')
    end
    
    ThetaP = repmat(theta, trajectoryCount, 1)';    
    X0P = repmat(x0, trajectoryCount, 1)';
else
    ThetaP = theta;
    X0P = x0;
end

    switch nargout
        case 0
            warning('No output specified')
            return
        case 1
            [Trajectories] = fTauleaping(X0P, ThetaP, modelStopTime, timeIntervals, opts);
        case 2
            [Trajectories, timeOut] = fTauleaping(X0P, ThetaP, modelStopTime, timeIntervals, opts);
        case 3
            [Trajectories, timeOut, R] = fTauleaping(X0P, ThetaP, modelStopTime, timeIntervals, opts);
        case 4
            [Trajectories, timeOut, R, IntG] = fTauleaping(X0P, ThetaP, modelStopTime, timeIntervals, opts);
        case 5
            [Trajectories, timeOut, R, IntG, simStats] = fTauleaping(X0P, ThetaP, modelStopTime, timeIntervals, opts);
        otherwise
            error('Too many output arguments specified')
    end
    Trajectories = permute(Trajectories, [2,1,3]);
end

