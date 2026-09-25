function sampledPredGlob = samplePredictedGlobal(globTable, globLogW, trackBernPred, detStates, hasState, numTracks)
%SAMPLEPREDICTEDGLOBAL Sample a predicted hypothesis given detected states.
% Constant columns cancel; score each used local hypothesis once per track.

if isempty(globLogW) || numTracks == 0
    sampledPredGlob = [];
    return
end

logPost = globLogW(:);

% Shared additive constants from constant columns do not affect sampling.
actIdx = getNonConstantTrackIdx(globTable, numTracks);
if isempty(actIdx)
    % All tracks constant across global hypotheses; only globLogW matters.
    aSample = sampleLogCategorical(logPost);
    sampledPredGlob = globTable(aSample, 1:numTracks);
    return
end

% Accumulate per-track contributions using a local-hypothesis lookup.
for kk = 1:numel(actIdx)
    i = actIdx(kk);

    col = globTable(:, i);
    idx = col + 1;  % hi=0 -> 1

    % Get unique hypothesis indices appearing in this column
    % (including possibly 0).
    u = unique(col);

    Hi = numel(trackBernPred{i});
    Li = -inf(Hi + 1, 1);

    if ~hasState(i)
        % hi==0 contributes nothing
        Li(1) = 0;

        % Only compute for unique positive hi values
        uPos = u(u > 0);
        for t = 1:numel(uPos)
            hi = uPos(t);
            lh = trackBernPred{i}{hi};
            % missed-detection term for this local hypothesis
            Li(hi + 1) = log1p(-lh.existence);
        end

    else
        % hasState(i)==true: hi==0 is invalid => -inf
        Li(1) = -inf;
        obj = detStates{i};

        % Only compute for unique positive hi values
        uPos = u(u > 0);
        for t = 1:numel(uPos)
            hi = uPos(t);
            lh = trackBernPred{i}{hi};
            Li(hi + 1) = detectedContribution(lh, obj);
        end
    end

    % Vectorized add over all global hypotheses
    logPost = logPost + Li(idx);
end

aSample = sampleLogCategorical(logPost);
sampledPredGlob = globTable(aSample, 1:numTracks);

    function acc = detectedContribution(lh, obj)
        % Full detected-state log contribution:
        % log(r) + state likelihood under GGIW.

        g = lh.ggiw;

        % Start with log existence
        r = lh.existence;
        if r <= 0
            acc = -inf;
            return
        end
        acc = log(r) + ggiwStateLogLik(g, obj);
    end

end
