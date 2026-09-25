function clutteredMeasurements = generateClutteredMeasurements(objectTracks, objectExtents, objectRates, parameters)
%GENERATECLUTTEREDMEASUREMENTS Generate Poisson object measurements and uniform clutter.
meanClutter = parameters.meanClutter;
region      = parameters.region;

numSteps   = size(objectTracks, 2);

xMin = region(1,1); xMax = region(2,1);
yMin = region(1,2); yMax = region(2,2);
xW = xMax - xMin;
yW = yMax - yMin;

clutteredMeasurements = cell(numSteps, 1);

for step = 1:numSteps
    validMask = ~isnan(objectTracks(1, step, :));
    objIdx = find(validMask(:));
    K = numel(objIdx);

    if K > 0
        counts = poissrnd(objectRates(step, objIdx));
        counts = counts(:);
    else
        counts = zeros(0,1);
    end
    totalObjMeas = sum(counts);

    numClutter = poissrnd(meanClutter);

    totalMeas = totalObjMeas + numClutter;
    meas = zeros(2, totalMeas);

    p = 1;
    for k = 1:K
        c = counts(k);
        if c == 0, continue; end

        o = objIdx(k);
        mu = objectTracks(1:2, step, o);
        C  = objectExtents(:, :, step, o);

        L = chol(C + 1e-9*eye(size(C)), 'lower');
        meas(:, p:(p+c-1)) = mu + L * randn(2, c);
        p = p + c;
    end

    if numClutter > 0
        q1 = totalObjMeas + 1;
        q2 = totalMeas;
        meas(1, q1:q2) = xW * rand(1, numClutter) + xMin;
        meas(2, q1:q2) = yW * rand(1, numClutter) + yMin;
    end

    if totalMeas > 1
        meas = meas(:, randperm(totalMeas));
    end

    clutteredMeasurements{step} = meas;
end
end
