function data = generateScenario(cfg, par, seed)
%GENERATESCENARIO Generate truth and measurements using separately seeded streams.
if isfield(par, 'truthSeed') && ~isempty(par.truthSeed)
    rng(par.truthSeed)
else
    rng default
end

[startStates, startExtents, startRates] = ...
    getStartStates(cfg.numObjects, cfg.startRadius, cfg.startSpeed, par);

[truthTracks, truthExtents, truthRates] = ...
    generateTracksUnknown(par, startStates, startExtents, startRates, cfg.appearance, cfg.numSteps);

if nargin < 3 || isempty(seed)
    seed = 0;
end
rng(seed)
measCell = generateClutteredMeasurements(truthTracks, truthExtents, truthRates, par);

data = struct();
data.startStates  = startStates;
data.startExtents = startExtents;
data.startRates   = startRates;

data.truthTracks  = truthTracks;
data.truthExtents = truthExtents;
data.truthRates   = truthRates;

data.measCell     = measCell;
data.appearance   = cfg.appearance;
data.numObjects   = cfg.numObjects;
data.numSteps     = cfg.numSteps;
end
