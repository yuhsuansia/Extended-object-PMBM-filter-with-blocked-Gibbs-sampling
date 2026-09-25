function measAssoc = sampleDataAssociation(Z, detStates, hasState, numTracks, numMeas, numSlots, cache)
%SAMPLEDATAASSOCIATION Sample associations conditional on sampled object existence.

logLambda = -inf(numSlots,1);
Hx        = zeros(2,numSlots);
S11       = zeros(numSlots,1);
S12       = zeros(numSlots,1);
S22       = zeros(numSlots,1);

H = cache.H;

% Precompute per-slot quantities only where hasState is true
for s = 1:numSlots
    if ~hasState(s), continue; end
    obj = detStates{s};

    logLambda(s) = log(obj.poissonRate);
    Hx(:,s)      = H * obj.kinematicState;

    S = obj.extentState;
    S11(s) = S(1,1); S12(s) = S(1,2); S22(s) = S(2,2);
end

% Indices of all existing slots (tracks + PPP slots) once
existIdx = find(hasState(:));
existTrackIdx = existIdx(existIdx <= numTracks);
existPPPIdx   = existIdx(existIdx > numTracks);

measAssoc = zeros(numMeas,1,'uint32');
for j = 1:numMeas
    z = Z(:,j);
    selfSlot = numTracks + j;

    if ~hasState(selfSlot)
        cand = [existTrackIdx; existPPPIdx(existPPPIdx > selfSlot)];
        if ~isempty(cand)
            ll = loggauss2_many(z, Hx(:,cand), S11(cand), S12(cand), S22(cand));
            ids = [cand; selfSlot];
            logw = [logLambda(cand) + ll; cache.logClutter];
            measAssoc(j) = uint32(sampleLogCategoricalSubset(ids, logw));
        else
            measAssoc(j) = uint32(selfSlot);
        end
    else
        measAssoc(j) = uint32(selfSlot);
    end
end
end
