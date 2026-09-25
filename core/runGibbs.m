function [trackBernUpd, globTableUpd] = runGibbs( ...
    trackBernUpd, hypoCache, currentGlob, trackBernPred, poissonUndetected, ...
    globTable, globLogW, Z, numTracks, numMeas, numSlots, par, cache)
%RUNGIBBS Alternate state, predicted-hypothesis, and association sampling.
% The first stored hypothesis is the initializer; the remaining rows are
% obtained from numGibbsIter - 1 Gibbs sweeps.

globTableUpd = zeros(par.numGibbsIter, numSlots, 'uint32');
globTableUpd(1,:) = currentGlob;

for r = 2:par.numGibbsIter
    detStates = sampleDetectedStatesFromGlobal( ...
        trackBernUpd, currentGlob, par.collapsedGibbs);

    if par.collapsedGibbs
        sampledPredGlob = samplePredictedGlobalCollapsed( ...
            globTable, globLogW, trackBernPred, detStates, numTracks);
        measAssoc = sampleDataAssociationCollapsed( ...
            Z, detStates, numTracks, numMeas, numSlots, cache);
    else
        hasState = ~cellfun(@isempty, detStates);
        sampledPredGlob = samplePredictedGlobal( ...
            globTable, globLogW, trackBernPred, detStates, hasState, numTracks);
        measAssoc = sampleDataAssociation( ...
            Z, detStates, hasState, numTracks, numMeas, numSlots, cache);
    end

    [trackBernUpd, hypoCache, currentGlob] = rebuildGlobalFromAssoc( ...
        trackBernUpd, hypoCache, trackBernPred, sampledPredGlob, poissonUndetected, ...
        Z, measAssoc, numTracks, numMeas, numSlots, par);
    globTableUpd(r,:) = currentGlob;
end
end
