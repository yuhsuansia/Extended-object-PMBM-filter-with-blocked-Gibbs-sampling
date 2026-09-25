function results = runFilter(data, par, verbose)
%RUNFILTER Run the PMBM filter on one supplied measurement sequence.
% data.measCell{k} is a 2-by-M_k matrix; use zeros(2,0) for an empty scan.
% Ground truth is used only by the GOSPA/PGOSPA evaluation below.
% results.state is predicted to the time step after the final measurement.
if nargin < 3
    verbose = true;
end
validateattributes(par.numGibbsIter, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
numSteps = numel(data.measCell);
assert(numSteps > 0, 'runFilter:EmptySequence', ...
    'Provide at least one measurement scan.');

cache = buildCache(par);
birth = buildPoissonBirth(par);
state = initFilterState(birth);
est = cell(numSteps, 1);
tStep = zeros(numSteps, 1);
pgospa = zeros(numSteps, 1);
pgospa_decomp(1, numSteps) = struct('localisation', 0, ...
    'existence_mismatch', 0, 'missed', 0, 'false', 0);
pG = 1; cG = 20; alphaG = 2;

if verbose
    fprintf('PMBM-Gibbs processing: ');
end
for k = 1:numSteps
    if verbose
        fprintf('%d ', k);
        if k == numSteps, fprintf('\n'); end
    end
    tk = tic;
    [state, est{k}, posterior] = pmbmStep(state, birth, data.measCell{k}, par, cache);
    tStep(k) = toc(tk);
    [pgospa(k), pgospa_decomp(k)] = ...
        computePGOSPAStepPMBM(data, k, posterior, pG, cG, alphaG);
end

[gospa, gospa_decomp] = computeGOSPAOverTime(data, est, pG, cG, alphaG);
results = struct('est', {est}, 'state', state, 'tStep', tStep, ...
    'totalFilterTime', sum(tStep), 'gospa', gospa, 'gospa_decomp', gospa_decomp, ...
    'pgospa', pgospa, 'pgospa_decomp', pgospa_decomp);
end
