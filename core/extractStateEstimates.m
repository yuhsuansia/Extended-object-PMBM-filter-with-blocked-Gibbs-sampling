function est = extractStateEstimates(trackBern, globTable, globLogW, existThreshold)
%EXTRACTSTATEESTIMATES Extract sufficiently probable objects from the MAP hypothesis.
if nargin < 4
    existThreshold = 0.5;
end

est = struct('slot',{},'existence',{}, ...
    'measurementRateMean',{},'mean',{},'extentMean',{});

if isempty(globLogW)
    return
end

[~, aMAP] = max(globLogW);
gMAP = globTable(aMAP, :);

K = numel(trackBern);
if K == 0
    return
end

numEst = 0;
% At most one estimate survives per slot under the MAP global hypothesis.
est(K) = struct('slot',[],'existence',[], ...
    'measurementRateMean',[],'mean',[],'extentMean',[]);

for s = 1:K
    hIdx = gMAP(s);
    if hIdx == 0
        continue
    end

    pool = trackBern{s};
    hypo = pool{hIdx};
    r    = hypo.existence;

    if r <= existThreshold
        continue
    end

    g = hypo.ggiw;
    extentMean = g.V / (g.v - 3);
    measurementRateMean = g.alpha / g.beta;

    numEst = numEst + 1;
    est(numEst) = struct( ...
        'slot',        s, ...
        'existence',   r, ...
        'measurementRateMean', measurementRateMean, ...
        'mean',        g.mean, ...
        'extentMean',  extentMean );
end

est = est(1:numEst);
end
