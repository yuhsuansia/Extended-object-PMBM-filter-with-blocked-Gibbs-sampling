function measAssoc = sampleDataAssociationCollapsed(Z, detStates, numTracks, numMeas, numSlots, cache)
%SAMPLEDATAASSOCIATIONCOLLAPSED Sample associations with marginalized existence probabilities.

H = cache.H;

logLambda = -inf(numSlots,1);
Hx        = nan(2,numSlots);
S11       = nan(numSlots,1);
S12       = nan(numSlots,1);
S22       = nan(numSlots,1);

slotExist = zeros(numSlots,1);

for s = 1:numSlots
    obj = detStates{s};
    r = obj.existence;
    slotExist(s) = r;

    if r > 0
        logLambda(s) = log(obj.poissonRate);
        Hx(:,s)      = H * obj.kinematicState;

        S = obj.extentState;
        S11(s) = S(1,1); S12(s) = S(1,2); S22(s) = S(2,2);
    end
end

existIdx = find(slotExist > 0); 
existTrackIdx = existIdx(existIdx <= numTracks);
existPPPIdx   = existIdx(existIdx > numTracks);

measAssoc = zeros(numMeas,1,'uint32');
for j = 1:numMeas
    z = Z(:,j);
    selfSlot = numTracks + j;
    rSelf = slotExist(selfSlot);

    if rSelf < 1
        cand = [existTrackIdx; existPPPIdx(existPPPIdx > selfSlot)];

        if ~isempty(cand)
            ll = loggauss2_many(z, Hx(:,cand), S11(cand), S12(cand), S22(cand));
            logwCand = logLambda(cand) + ll + log1p(-rSelf) + log(slotExist(cand));
            ids = [cand; selfSlot];
            if rSelf <= 0
                logwSelf = cache.logClutter;
            else
                llSelf = loggauss2_many(z, Hx(:,selfSlot), S11(selfSlot), S12(selfSlot), S22(selfSlot));
                a = logLambda(selfSlot) + llSelf + log(rSelf);
                b = cache.logClutter + log1p(-rSelf);
                logwSelf = logsumexp2(a, b);
            end
            measAssoc(j) = uint32(sampleLogCategoricalSubset(ids, [logwCand; logwSelf]));
        else
            measAssoc(j) = uint32(selfSlot);
        end
    else
        measAssoc(j) = uint32(selfSlot);
    end
end

measAssoc = reorderNewBernoulliLabels(measAssoc, numTracks, numMeas, numSlots);

end
