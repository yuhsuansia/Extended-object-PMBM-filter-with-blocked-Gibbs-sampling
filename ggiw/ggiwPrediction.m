function predictedGGIW = ggiwPrediction(updatedGGIW, parameters, cache)
%GGIWPREDICTION Predict kinematics, gamma rate, and inverse-Wishart extent.
dt   = parameters.dt;
eta  = parameters.gammaEta;
tau  = parameters.iwTau;

ff = exp(-dt / tau);

A = cache.A;
Q = cache.Q;

a = updatedGGIW.alpha;
b = updatedGGIW.beta;
m = updatedGGIW.mean;
P = updatedGGIW.cov;
v = updatedGGIW.v;
V = updatedGGIW.V;

predictedGGIW.alpha = a / eta;
predictedGGIW.beta  = b / eta;

predictedGGIW.mean = A * m;
predictedGGIW.cov  = A * P * A.' + Q;

predictedGGIW.v = 3 + ff * (v - 3);
predictedGGIW.V = ff * V;

predictedGGIW = ggiwFinalizeCache(predictedGGIW);
end
