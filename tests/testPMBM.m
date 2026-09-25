function tests = testPMBM
%TESTPMBM Deterministic checks for the standalone PMBM implementation.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
testCase.TestData.originalRng = rng;
rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
setupPath;
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
rng(testCase.TestData.originalRng);
end

function testScenarioReproducibility(testCase)
cfg = defaultConfig();
cfg.numSteps = 16;
par = buildParameters(cfg);
first = generateScenario(cfg, par, cfg.measurementSeed);
cfg.numGibbsIter = 3;
cfg.collapsedGibbs = ~cfg.collapsedGibbs;
cfg.enableClustering = ~cfg.enableClustering;
second = generateScenario(cfg, buildParameters(cfg), cfg.measurementSeed);
verifyTrue(testCase, isequaln(first, second));
third = generateScenario(cfg, par, cfg.measurementSeed + 1);
verifyTrue(testCase, isequaln(first.truthTracks, third.truthTracks));
verifyFalse(testCase, isequal(first.measCell, third.measCell));
end

function testFilterModes(testCase)
cfg = defaultConfig();
cfg.numSteps = 12;
cfg.numGibbsIter = 5;
cfg.plotGroundTruth = false;
data = generateScenario(cfg, buildParameters(cfg), cfg.measurementSeed);
data.measCell{8} = zeros(2, 0);
data.measCell{9} = data.measCell{9}(:, 1);
for collapsed = [false, true]
    for clustered = [false, true]
        cfg.collapsedGibbs = collapsed;
        cfg.enableClustering = clustered;
        par = buildParameters(cfg);
        rng(123, 'twister');
        result = runFilter(data, par, false);
        verifySize(testCase, result.gospa, [cfg.numSteps, 1]);
        verifyTrue(testCase, all(isfinite(result.gospa)));
        verifyTrue(testCase, all(isfinite(result.pgospa)));
        components = [result.gospa_decomp.localisation] ...
            + [result.gospa_decomp.missed] + [result.gospa_decomp.false];
        verifyEqual(testCase, result.gospa, components.', 'AbsTol', 1e-10);
        verifyEqual(testCase, sum(exp(result.state.globLogW)), 1, 'AbsTol', 1e-12);
        verifyLessThanOrEqual(testCase, size(result.state.globTable, 1), cfg.numGibbsIter);
        rng(123, 'twister');
        repeated = runFilter(data, par, false);
        verifyTrue(testCase, isequaln(result.state, repeated.state));
        verifyEqual(testCase, result.est, repeated.est);
    end
end
end

function testEmptyScans(testCase)
cfg = defaultConfig();
cfg.numGibbsIter = 1;
par = buildParameters(cfg);
cache = buildCache(par);
birth = buildPoissonBirth(par);
state = initFilterState(birth);
for k = 1:3
    [state, estimates] = pmbmStep(state, birth, zeros(2, 0), par, cache);
    verifyEmpty(testCase, estimates);
    verifyEmpty(testCase, state.trackBern);
    verifyTrue(testCase, all(isfinite([state.poissonUndetected.logWeight])));
end
end

function testDbscanLabels(testCase)
cfg = defaultConfig();
par = buildParameters(cfg);
Z = [0, 0.1, 100; 0, 0.1, 100];
labels = dbscanAssignUnassigned(Z, zeros(3,1,'uint32'), 4, par);
verifyEqual(testCase, labels, uint32([6; 6; 7]));
end

function testGGIWUpdates(testCase)
par = buildParameters(defaultConfig());
birth = buildPoissonBirth(par);
g = birth.ggiw;
Z = [1, 2, -1; 0, -2, 1];
[updated, ll] = ggiwUpdate(g, Z, par);
verifyEqual(testCase, updated.alpha, g.alpha + 3);
verifyEqual(testCase, updated.beta, g.beta + 1);
verifyEqual(testCase, updated.v, g.v + 3);
verifyTrue(testCase, isfinite(ll));
verifyGreaterThan(testCase, min(eig(updated.cov)), 0);
verifyGreaterThan(testCase, min(eig(updated.V)), 0);
[singleton, singletonLL] = ggiwUpdate1(g, Z(:,1), par);
[general, generalLL] = ggiwUpdate(g, Z(:,1), par);
verifyEqual(testCase, singleton.mean, general.mean, 'AbsTol', 1e-12);
verifyEqual(testCase, singleton.V, general.V, 'AbsTol', 1e-12);
verifyEqual(testCase, singletonLL, generalLL, 'AbsTol', 1e-8);
newObject = pppUpdate(birth, Z, par);
verifyEqual(testCase, newObject.existence, 1);
singletonObject = pppUpdate(birth, Z(:,1), par);
verifyGreaterThan(testCase, singletonObject.existence, 0);
verifyLessThan(testCase, singletonObject.existence, 1);
end

function testCategoricalEdges(testCase)
verifyEqual(testCase, logsumexp([]), -inf);
verifyEqual(testCase, logsumexp([-inf; -inf]), -inf);
verifyError(testCase, @() sampleLogCategorical([-inf; -inf]), ...
    'sampleLogCategorical:AllImpossible');
verifyError(testCase, @() sampleLogCategorical(NaN), 'sampleLogCategorical:NaNWeights');
verifyEqual(testCase, sampleLogCategorical([-inf; 0; -inf]), 2);
before = rng;
verifyEqual(testCase, sampleLogCategoricalSubset(7, 0), 7);
verifyEqual(testCase, rng, before);
end
