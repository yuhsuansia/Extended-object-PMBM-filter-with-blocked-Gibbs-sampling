function [objectTracks, objectExtents, objectRates] = generateTracksUnknown(parameters, startStates, extentMatrixes, startMeasurementRates, appearanceFromTo, numSteps)
%GENERATETRACKSUNKNOWN Generate trajectories with specified birth and death times.
accDev = parameters.sigmaAcc;
dt     = parameters.dt;

numObjects = size(startStates, 2);

motionModel = getConstantVelocityModel(dt, accDev);
A = motionModel.A;
Q = motionModel.Q;

LQ = chol(Q, 'lower');

objectTracks  = nan(4, numSteps, numObjects);
objectExtents = nan(2, 2, numSteps, numObjects);
objectRates   = nan(numSteps, numObjects);

for obj = 1:numObjects
    x = startStates(:, obj);

    a0 = max(1, appearanceFromTo(1, obj));
    a1 = min(numSteps, appearanceFromTo(2, obj));

    if a0 > a1
        for t = 1:numSteps
            x = A*x + LQ*randn(4,1);
        end
        continue
    end

    for t = 1:(a0-1)
        x = A*x + LQ*randn(4,1);
    end

    for t = a0:a1
        x = A*x + LQ*randn(4,1);
        objectTracks(:, t, obj) = x;
    end

    for t = (a1+1):numSteps
        x = A*x + LQ*randn(4,1);
    end

    T = a1 - a0 + 1;
    objectExtents(:, :, a0:a1, obj) = repmat(extentMatrixes(:, :, obj), 1, 1, T);
    objectRates(a0:a1, obj) = startMeasurementRates(obj);
end
end
