function [trackBernUpd, hypoCache, currentGlob] = rebuildGlobalFromAssoc( ...
    trackBernUpd, hypoCache, trackBernPred, sampledPredGlob, poissonUndetected, ...
    Z, measAssoc, numTracks, numMeas, numSlots, par)
%REBUILDGLOBALFROMASSOC Rebuild local/global hypotheses from sampled association.

groups = groupMeasurementsByAssoc(measAssoc, numSlots);

currentGlob = zeros(1, numSlots, 'uint32');

% existing tracks
for i = 1:numTracks
    parentIdx = sampledPredGlob(i);
    if parentIdx == 0, continue; end

    idxMeas  = groups{i};
    measIdx  = uint16(idxMeas(:).');
    key      = makeHypoKey(parentIdx, measIdx, true);

    % cache check BEFORE constructing upd ---
    [hIdx, found] = getLocalHypoIdx(hypoCache{i}, key);
    if found
        currentGlob(i) = hIdx;
        continue
    end

    % build only if needed
    predHypo = trackBernPred{i}{parentIdx};
    if isempty(idxMeas)
        upd = bernoulliMissedDetection(predHypo);
    else
        upd = bernoulliUpdate(predHypo, Z(:,idxMeas), par);
    end
    upd.parentalIdx = parentIdx;
    upd.measIdx     = measIdx;

    [trackBernUpd{i}, hypoCache{i}, hIdx] = addLocalHypo(trackBernUpd{i}, hypoCache{i}, key, upd);
    currentGlob(i) = hIdx;
end

% PPP slots
for j = 1:numMeas
    s = numTracks + j;
    idxMeas = groups{s};

    if isempty(idxMeas)
        currentGlob(s) = uint32(0);
        continue
    end

    measIdx = uint16(idxMeas(:).');
    key     = makeHypoKey(0, measIdx, true);

    % cache check BEFORE constructing upd ---
    [hIdx, found] = getLocalHypoIdx(hypoCache{s}, key);
    if found
        currentGlob(s) = hIdx;
        continue
    end

    % build only if needed
    upd = pppUpdate(poissonUndetected, Z(:,idxMeas), par);
    upd.parentalIdx = 0;
    upd.measIdx     = measIdx;

    [trackBernUpd{s}, hypoCache{s}, hIdx] = addLocalHypo(trackBernUpd{s}, hypoCache{s}, key, upd);
    currentGlob(s) = hIdx;
end
end
