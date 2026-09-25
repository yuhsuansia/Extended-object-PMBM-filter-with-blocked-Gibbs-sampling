function measAssoc = greedyInitAssocMeanSurrogate( ...
    Z, numMeas, numTracks, bestGlob, trackBern, poissonUndetected, par)
%GREEDYINITASSOCMEANSURROGATE Greedy initializer using mean-based GGIW scores.
%
% Existing-track score:
%   log r + log E[lambda] + log N(z | H x, H P H' + E[X])
%
% PPP score:
%   log intensity + log E[lambda] + log N(z | H x, H P H' + E[X])
% marginalized with clutter.

measAssoc = zeros(numMeas, 1, 'uint32');
H = par.H;
logClutter = log(par.clutterIntensity);

for j = 1:numMeas
    z = Z(:, j);
    greedyLL = -inf(numTracks + 1, 1);

    for i = 1:numTracks
        parentIdx = bestGlob(i);
        if parentIdx == 0
            continue
        end

        predHypo = trackBern{i}{parentIdx};
        r = predHypo.existence;
        if r <= 0
            continue
        end

        g = predHypo.ggiw;
        lambdaMean = g.alpha / g.beta;
        ll = meanSurrogateLogLik(g, z, H);
        greedyLL(i) = log(r) + log(lambdaMean) + ll;
    end

    numPPP = numel(poissonUndetected);
    logWeights = -inf(numPPP, 1);
    for ii = 1:numPPP
        g = poissonUndetected(ii).ggiw;
        lambdaMean = g.alpha / g.beta;
        ll = meanSurrogateLogLik(g, z, H);
        logWeights(ii) = poissonUndetected(ii).logWeight + log(lambdaMean) + ll;
    end

    greedyLL(end) = logsumexp2(logsumexp(logWeights), logClutter);

    [~, best] = max(greedyLL);
    if best <= numTracks
        measAssoc(j) = uint32(best);
    else
        % DBSCAN groups the measurements left unassigned by this pass.
        if par.enableClustering
            measAssoc(j) = uint32(0);
        else
            measAssoc(j) = uint32(numTracks + j);
        end
    end
end
end

function ll = meanSurrogateLogLik(g, z, H)
% Approximate the predictive measurement density with the means of the
% gamma and inverse-Wishart factors, while keeping the kinematic
% uncertainty contribution H P H'.
Xmean = getExtentMean(g);
HP = H * g.cov;
S = HP * H.' + Xmean;
S = 0.5 * (S + S.');
mu = H * g.mean;
ll = loggauss2_many(z, mu, S(1,1), S(1,2), S(2,2));
end

function Xmean = getExtentMean(g)
% Use cached E[X] when available, otherwise fall back to V / (v - 3).
if isfield(g, 'Xhat') && ~isempty(g.Xhat)
    Xmean = g.Xhat;
else
    Xmean = g.V / (g.v - 3);
    Xmean = 0.5 * (Xmean + Xmean.');
end
end
