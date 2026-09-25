function sampledPredGlob = samplePredictedGlobalCollapsed(globTable, globLogW, trackBernPred, detStates, numTracks)
%SAMPLEPREDICTEDGLOBALCOLLAPSED Sample a predicted hypothesis, marginalizing existence.
% Constant columns cancel; score each used local hypothesis once per track.

if isempty(globLogW) || numTracks == 0
    sampledPredGlob = [];
    return
end

logPost = globLogW(:);

% Precompute detected existence and log(1-rdet).
rDet = zeros(numTracks, 1);
for i = 1:numTracks
    rDet(i) = detStates{i}.existence;
end

log1m_rDet = log1p(-rDet);

% Shared additive constants from constant columns do not affect sampling.
actIdx = getNonConstantTrackIdx(globTable, numTracks);
if isempty(actIdx)
    % All tracks constant across global hypotheses: only globLogW matters.
    aSample = sampleLogCategorical(logPost);
    sampledPredGlob = globTable(aSample, 1:numTracks);
    return
end

% Accumulate per-track contributions using a local-hypothesis lookup.
for kk = 1:numel(actIdx)
    i = actIdx(kk);

    col = globTable(:, i);
    idx = col + 1;            % hi=0 -> 1

    u = unique(col);          % unique hi appearing for this track
    Hi = numel(trackBernPred{i});

    % Lookup for hi=0..Hi. Initialize to -inf; we'll fill the ones we need.
    Li = -inf(Hi + 1, 1);

    % A missing predicted track contributes its detected nonexistence term.
    Li(1) = log1m_rDet(i);

    % With no detected existence, only predicted nonexistence contributes.
    rdet = rDet(i);
    if rdet <= 0
        uPos = u(u > 0);
        for t = 1:numel(uPos)
            hi = uPos(t);
            r  = trackBernPred{i}{hi}.existence;
            Li(hi + 1) = log1p(-r);
        end

        % Vectorized add
        logPost = logPost + Li(idx);
        continue
    end

    % Otherwise: rdet > 0, for hi>0 we need the full collapsed existence term:
    % termExist = logsumexp( log(rdet)+log(r)+logLik , log(1-rdet)+log(1-r) )
    obj = detStates{i};
    log_rdet = log(rdet);
    log1m_rdet = log1m_rDet(i);

    uPos = u(u > 0);
    for t = 1:numel(uPos)
        hi = uPos(t);
        lh = trackBernPred{i}{hi};
        r  = lh.existence;

        if r <= 0
            % log(1-r) ~ 0 when r=0; and first branch is -inf
            Li(hi + 1) = log1m_rdet + log1p(-r);
            continue
        end

        logLik = ggiwStateLogLik(lh.ggiw, obj);

        % two branches
        a1 = log_rdet + log(r) + logLik;
        a2 = log1m_rdet + log1p(-r);

        Li(hi + 1) = logsumexp2(a1, a2);
    end

    % Vectorized add over all global hypotheses
    logPost = logPost + Li(idx);
end

aSample = sampleLogCategorical(logPost);
sampledPredGlob = globTable(aSample, 1:numTracks);

end
