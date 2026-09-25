function par = buildParameters(cfg)
%BUILDPARAMETERS Model, initialization, and hypothesis-pruning parameters.
validateattributes(cfg.numGibbsIter, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
validateattributes(cfg.meanMeasPerObject, {'numeric'}, ...
    {'scalar', 'positive', 'finite'});
validateattributes(cfg.meanClutter, {'numeric'}, ...
    {'scalar', 'positive', 'finite'});
par = struct();

par.dt         = 0.2;
par.sigmaAcc   = 0.8;
par.gammaEta   = 1.01;
par.iwTau      = par.dt * 100;
par.pSurvive   = 0.99;
par.meanBirths = 0.01;

par.region = [-150 -150; 150 150];
regionSize = par.region(2,:) - par.region(1,:);

par.meanMeasPerObject = cfg.meanMeasPerObject;
par.H = [1 0 0 0; 0 1 0 0];

par.meanClutter      = cfg.meanClutter;
par.clutterIntensity = par.meanClutter / prod(regionSize);
par.truthSeed        = cfg.truthSeed;

par.extentNominal = cfg.extentNominal;
par.iwDoF = 100;
par.iwV0  = cfg.extentNominal * (par.iwDoF - 3);
par.birthIwDoF = 100;

par.gammaScale = 100;
par.gammaShape = par.meanMeasPerObject * par.gammaScale;
par.birthGammaScale = 100;

par.enableClustering = cfg.enableClustering;
par.dbscanEps    = sqrt(trace(par.extentNominal)/2);
par.dbscanMinPts = 2;

par.numGibbsIter                 = cfg.numGibbsIter;
par.globalHypoPruneThreshold     = 1e-3;
par.localHypoExistenceThreshold  = 1e-3;
par.poissonUndetectedPruneThresh = 1e-3;
par.existThreshold               = 0.5;

par.collapsedGibbs = cfg.collapsedGibbs;

end
