function [updatedGGIW, predictedLogLik] = ggiwUpdate(predictedGGIW, measurements, parameters)
%GGIWUPDATE Update a GGIW density with a nonempty measurement cell.
H = parameters.H;
m = size(measurements, 2);

z_bar = sum(measurements, 2) / m;
% Within-group scatter used by the multi-measurement extent update.
Z = measurements*measurements.' - m*(z_bar*z_bar.');
Z = 0.5*(Z+Z.');

updatedGGIW.alpha = predictedGGIW.alpha + m;
updatedGGIW.beta  = predictedGGIW.beta  + 1;

% Reuse the cached IW mean extent instead of rebuilding V/(v-3) each call.
Xhat = predictedGGIW.Xhat;

P  = predictedGGIW.cov;
HP = H * P;
S  = HP * H.' + Xhat / m;

S = 0.5*(S+S.');
LS = chol(S + 1e-9*eye(size(S)), 'lower');

Y = LS \ HP;
Y = LS.' \ Y;
K = Y.';

eps = z_bar - H * predictedGGIW.mean;

updatedGGIW.mean = predictedGGIW.mean + K * eps;

Pnew = P - K * HP;
updatedGGIW.cov = 0.5 * (Pnew + Pnew.');

u = LS \ eps;
u = LS.' \ u;

% N captures the centroid innovation, while Z captures within-group scatter.
Xu = Xhat * u;
N  = Xu * Xu.';

updatedGGIW.v = predictedGGIW.v + m;
Vnew = predictedGGIW.V + N + Z;
updatedGGIW.V = 0.5 * (Vnew + Vnew.');

logdetS = 2 * sum(log(diag(LS)));

LVupd = chol(updatedGGIW.V, 'lower');
logdetVupd = 2 * sum(log(diag(LVupd)));

predictedLogLik = -m * log(pi) - log(m);
predictedLogLik = predictedLogLik ...
    + (predictedGGIW.v/2) * predictedGGIW.logdetV ...
    - (updatedGGIW.v/2)   * logdetVupd;
predictedLogLik = predictedLogLik ...
    + gammaln2(updatedGGIW.v/2) - gammaln2(predictedGGIW.v/2);
predictedLogLik = predictedLogLik + 0.5*(predictedGGIW.logdetXhat - logdetS);
predictedLogLik = predictedLogLik ...
    + gammaln(updatedGGIW.alpha) - gammaln(predictedGGIW.alpha);
predictedLogLik = predictedLogLik ...
    + predictedGGIW.alpha * log(predictedGGIW.beta) ...
    - updatedGGIW.alpha   * log(updatedGGIW.beta);

updatedGGIW = ggiwFinalizeCache(updatedGGIW);
end
