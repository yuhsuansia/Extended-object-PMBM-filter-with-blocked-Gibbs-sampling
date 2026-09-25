function [startStates, startMatrixes, startMeasurementRates] = getStartStates(numObjects, radius, speed, parameters)
%GETSTARTSTATES Sample initial extents/rates and arrange inward-moving objects.
nu  = parameters.iwDoF;
Psi = parameters.iwV0;

startMatrixes = zeros(2,2,numObjects);
for k = 1:numObjects
    startMatrixes(:,:,k) = iwishrnd(Psi, nu);
end

startMeasurementRates = gamrnd(parameters.gammaShape, 1/parameters.gammaScale, numObjects, 1);

if numObjects < 2
    startStates = zeros(4,1);
    startStates(3) = speed;
    return
end

angles = (0:numObjects-1) * (2*pi/numObjects);
sa = sin(angles);
ca = cos(angles);

startStates = [ ...
    sa * radius; ...
    ca * radius; ...
    -sa * speed; ...
    -ca * speed ];

startStates(:,1) = [0; radius; 0; -speed];
end
