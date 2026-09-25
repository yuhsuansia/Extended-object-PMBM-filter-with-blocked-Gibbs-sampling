function [updatedGGIW, logLik] = ggiwUpdate1(predGGIW, z, parameters)
%GGIWUPDATE1 Singleton PPP update with its dedicated numerical regularization.
H = parameters.H;

alpha = predGGIW.alpha;
beta  = predGGIW.beta;
v     = predGGIW.v;
V     = predGGIW.V;
m0    = predGGIW.mean;
P0    = predGGIW.cov;

P0 = 0.5*(P0 + P0.');
V  = 0.5*(V  + V.');

% Reuse cached IW mean extent and log-determinants from ggiwFinalizeCache.
Xhat = predGGIW.Xhat;

eps = z - H*m0;

HP = H*P0;
S  = HP*H.' + Xhat;
S  = 0.5*(S + S.');

LS = chol(S + 1e-9*eye(2), 'lower');
logdetS = 2*sum(log(diag(LS)));

Y = LS \ HP;
Y = LS.' \ Y;
K = Y.';

m1 = m0 + K*eps;
P1 = P0 - K*HP;
P1 = 0.5*(P1 + P1.');

u = LS \ eps;
u = LS.' \ u;

% Rank-1 extent innovation from the centered measurement residual.
Xu = Xhat * u;
N  = Xu * Xu.';
N  = 0.5*(N + N.');

v1 = v + 1;
V1 = V + N;
V1 = 0.5*(V1 + V1.');

alpha1 = alpha + 1;
beta1  = beta  + 1;

logLik = -log(pi);

LV1 = chol(V1 + 1e-9*eye(2), 'lower');
logdetV1 = 2*sum(log(diag(LV1)));

logLik = logLik ...
    + (v/2)   * predGGIW.logdetV ...
    - (v1/2)  * logdetV1 ...
    + gammaln2(v1/2) - gammaln2(v/2) ...
    + 0.5*(predGGIW.logdetXhat - logdetS);

logLik = logLik ...
    + gammaln(alpha1) - gammaln(alpha) ...
    + alpha*log(beta) - alpha1*log(beta1);

updatedGGIW = predGGIW;
updatedGGIW.alpha = alpha1;
updatedGGIW.beta  = beta1;
updatedGGIW.mean  = m1;
updatedGGIW.cov   = P1;
updatedGGIW.v     = v1;
updatedGGIW.V     = V1;

updatedGGIW = ggiwFinalizeCache(updatedGGIW);
end
