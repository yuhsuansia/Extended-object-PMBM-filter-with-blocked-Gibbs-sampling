%MAIN Run one extended-object PMBM-Gibbs filtering demonstration.
clearvars; close all; rng default
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setupPath;

% Edit defaultConfig for the scenario, sampler, initialization, and seeds.
cfg = defaultConfig();
par = buildParameters(cfg);
data = generateScenario(cfg, par, cfg.measurementSeed);

if cfg.plotGroundTruth
    plotGroundTruth(data, par);
end

results = runFilter(data, par);

fprintf('\nTiming:\n');
fprintf('  total filtering time = %.3f s\n', results.totalFilterTime);
fprintf('  mean per step        = %.3f ms\n', 1e3 * mean(results.tStep));
fprintf('\nSummary:\n');
fprintf('  mean(GOSPA) = %.4g\n', mean(results.gospa));
fprintf('  mean(PGOSPA) = %.4g\n', mean(results.pgospa));
fprintf('  mean(Localisation) = %.4g\n', mean([results.gospa_decomp.localisation]));
fprintf('  mean(Missed)       = %.4g\n', mean([results.gospa_decomp.missed]));
fprintf('  mean(False)        = %.4g\n', mean([results.gospa_decomp.false]));
fprintf('  mean(PG Localisation) = %.4g\n', mean([results.pgospa_decomp.localisation]));
fprintf('  mean(PG ExistenceMismatch) = %.4g\n', mean([results.pgospa_decomp.existence_mismatch]));
fprintf('  mean(PG Missed)       = %.4g\n', mean([results.pgospa_decomp.missed]));
fprintf('  mean(PG False)        = %.4g\n', mean([results.pgospa_decomp.false]));
