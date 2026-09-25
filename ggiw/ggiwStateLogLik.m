function logLik = ggiwStateLogLik(g, obj)
%GGIWSTATELOGLIK Log-likelihood of sampled state under a cached GGIW model.

lambda = obj.poissonRate;
if lambda <= 0
    logLik = -inf;
    return
end
logLik = g.gammaConst + g.alphaMinus1 * log(lambda) - g.beta * lambda;

dx = obj.kinematicState - g.mean;
y  = g.covL \ dx;
logLik = logLik + g.logNConst - 0.5 * (y.' * y);

X = obj.extentState;
x11 = X(1,1); x12 = X(1,2); x22 = X(2,2);

detX = x11*x22 - x12*x12;
if detX <= 0
    detX = realmin;
end
logdetX = log(detX);

inv11 =  x22 / detX;
inv22 =  x11 / detX;
inv12 = -x12 / detX;

Psi = g.V;
trTerm = Psi(1,1)*inv11 + 2*Psi(1,2)*inv12 + Psi(2,2)*inv22;

logLik = logLik + g.iwConst - g.iwCoef * logdetX - 0.5 * trTerm;
end
