function [state, estk, postk] = pmbmStep(state, birth, Z, par, cache)
%PMBMSTEP One time-step: update (Gibbs) + prune + predict.

numMeas   = size(Z,2);
numTracks = numel(state.trackBern);
numSlots  = numTracks + numMeas;

% 1) Select initial predicted global hypothesis
bestGlob = selectBestGlobalHypo(state.globTable, state.globLogW, numTracks);

% 2) Prepare updated pools + caches
[trackBernUpd, hypoCache] = initUpdatedPools(numSlots);

% 3) Initial association, optionally clustering unassigned measurements.
measAssoc = greedyInitAssocMeanSurrogate(Z, numMeas, numTracks, bestGlob, state.trackBern, state.poissonUndetected, par);
if par.enableClustering
    measAssoc = dbscanAssignUnassigned(Z, measAssoc, numTracks, par);
end

% 4) Use the same hypothesis construction for initialization and sampling.
[trackBernUpd, hypoCache, currentGlob] = rebuildGlobalFromAssoc( ...
    trackBernUpd, hypoCache, state.trackBern, bestGlob, state.poissonUndetected, ...
    Z, measAssoc, numTracks, numMeas, numSlots, par);

% 5) Gibbs sampling.
[trackBernUpd, globTableUpd] = runGibbs( ...
    trackBernUpd, hypoCache, currentGlob, state.trackBern, state.poissonUndetected, ...
    state.globTable, state.globLogW, Z, numTracks, numMeas, numSlots, par, cache);

% 6) Deduplicate, weight, and prune hypotheses.
[trackBernUpd, globTableUpd, globLogWUpd] = postGibbsAndPrune(trackBernUpd, globTableUpd, par);

% Commit posterior
state.trackBern = trackBernUpd;
state.globTable = globTableUpd;
state.globLogW  = globLogWUpd;

% Estimates
estk = extractStateEstimates(state.trackBern, state.globTable, state.globLogW, par.existThreshold);

% Optional: posterior snapshot for downstream metrics (before prediction).
if nargout >= 3
    postk = struct();
    postk.trackBern = trackBernUpd;
    postk.globTable = globTableUpd;
    postk.globLogW  = globLogWUpd;
end

% 7) Predict the PPP and Bernoulli components to the next time step.
state.poissonUndetected = predictPPP(state.poissonUndetected, birth, par, cache);

state.trackBern = predictBernoullis(state.trackBern, par, cache);
end
