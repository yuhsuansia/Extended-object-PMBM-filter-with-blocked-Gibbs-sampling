function g = ggiwFinalizeCache(g)
% Cache all per-GGIW constants/factorizations needed in Gibbs scoring + sampling.

% --- Gaussian cache (4D) ---
g.cov = 0.5*(g.cov + g.cov.');
Lc = chol(g.cov + 1e-6*eye(4), 'lower');
g.covL = Lc;
g.logdetCov = 2*sum(log(diag(Lc)));
g.logNConst = -0.5*(4*log(2*pi) + g.logdetCov);

% --- Gamma cache ---
g.gammaConst   = g.alpha*log(g.beta) - gammaln(g.alpha);
g.alphaMinus1  = g.alpha - 1;

% --- IW cache (p=2) ---
g.V = 0.5*(g.V + g.V.');
Lv = chol(g.V + 1e-6*eye(2), 'lower');
g.logdetV = 2*sum(log(diag(Lv)));
% The update code repeatedly uses the IW mean extent and its log-determinant.
g.Xhat = g.V / (g.v - 3);
g.Xhat = 0.5*(g.Xhat + g.Xhat.');
g.logdetXhat = g.logdetV - 2*log(g.v - 3); % p=2
g.iwConst = (g.v/2)*g.logdetV - (g.v*2/2)*log(2) - gammaln2(g.v/2); % p=2
g.iwCoef  = (g.v + 3)/2;  % (v + p + 1)/2 with p=2

% cache chol(inv(V)) for fast IW sampling
invV = Lv'\(Lv\eye(2));
invV = 0.5*(invV+invV.');
g.VinvL = chol(invV + 1e-9*eye(2), 'lower');
end
