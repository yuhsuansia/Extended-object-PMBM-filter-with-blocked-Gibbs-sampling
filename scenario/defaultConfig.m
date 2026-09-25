function cfg = defaultConfig()
%DEFAULTCONFIG Settings for one reproducible PMBM-Gibbs demonstration.
cfg = struct();
cfg.numGibbsIter = 100;
cfg.collapsedGibbs = true;
cfg.enableClustering = false;
cfg.meanMeasPerObject = 10;
cfg.meanClutter = 10;
cfg.truthSeed = 0;
cfg.measurementSeed = 0;

cfg.numSteps      = 100;
cfg.numObjects    = 10;
cfg.extentNominal = diag([5 5]);
cfg.startRadius   = 125;
cfg.startSpeed    = 12.5;
cfg.plotGroundTruth = true;

% Birth/death times
cfg.appearance = ...
    [[3;83],[3;83], ...
    [6;86],[6;86], ...
    [9;89],[9;89], ...
    [12;92],[12;92], ...
    [15;95],[15;95]];
end
