function updatedLocalHypo = pppUpdate(predictedPPP, measurements, parameters)
%PPPUPDATE Form a new Bernoulli from a PPP mixture and a measurement cell.
numGGIW = numel(predictedPPP);
m = size(measurements, 2);
clutterI = parameters.clutterIntensity;
useSingletonUpdate = (m == 1);

if numGGIW == 1
    if useSingletonUpdate
        [ggiwUp, logLik] = ggiwUpdate1(predictedPPP.ggiw, measurements, parameters);
    else
        [ggiwUp, logLik] = ggiwUpdate(predictedPPP.ggiw, measurements, parameters);
    end

    logW = predictedPPP.logWeight + logLik;

    updatedLocalHypo.ggiw = ggiwUp;

    if useSingletonUpdate
        logDen = logsumexp2(logW, log(clutterI));
        updatedLocalHypo.existence = exp(logW - logDen);
        updatedLocalHypo.logWeight = logDen;
    else
        updatedLocalHypo.existence = 1;
        updatedLocalHypo.logWeight = logW;
    end
    return
end

alphas = zeros(numGGIW, 1);
betas  = zeros(numGGIW, 1);
means  = zeros(4, numGGIW);
covs   = zeros(4, 4, numGGIW);
vs     = zeros(numGGIW, 1);
Vs     = zeros(2, 2, numGGIW);
logWeights = zeros(numGGIW, 1);

for i = 1:numGGIW
    % Singleton groups can use the dedicated 1-measurement update cheaply.
    if useSingletonUpdate
        [gi, ll] = ggiwUpdate1(predictedPPP(i).ggiw, measurements, parameters);
    else
        [gi, ll] = ggiwUpdate(predictedPPP(i).ggiw, measurements, parameters);
    end

    alphas(i)      = gi.alpha;
    betas(i)       = gi.beta;
    means(:, i)    = gi.mean;
    covs(:, :, i)  = gi.cov;
    vs(i)          = gi.v;
    Vs(:, :, i)    = gi.V;

    logWeights(i) = predictedPPP(i).logWeight + ll;
end

logSumW = logsumexp(logWeights);

if useSingletonUpdate
    logDen = logsumexp2(logSumW, log(clutterI));
    updatedLocalHypo.existence = exp(logSumW - logDen);
    updatedLocalHypo.logWeight = logDen;
else
    updatedLocalHypo.existence = 1;
    updatedLocalHypo.logWeight = logSumW;
end

weights = exp(logWeights - logSumW);

ggiw = struct();
[ggiw.alpha, ggiw.beta] = gammaMerge(alphas, betas, weights);
[ggiw.mean,  ggiw.cov]  = gaussianMerge(means, covs, weights);
[ggiw.v,     ggiw.V]    = inverseWishartMerge(vs, Vs, weights);

ggiw = ggiwFinalizeCache(ggiw);

updatedLocalHypo.ggiw = ggiw;
end
