function [pgospa, pgospa_decomp] = computePGOSPAStepPMBM(data, k, postk, p, c, alpha)
%COMPUTEPGOSPASTEPPMBM PMBM weighted-sum PGOSPA at one time step.
%
% For each MB (global hypothesis), compute PGOSPA to the ground truth and
% weight it by the MB weight.

truth = buildTruthSet(data, k);
x_r = ones(1, size(truth.x, 2));

% Handle empty/degenerate posterior robustly as one empty MB.
if isempty(postk) || ~isfield(postk, 'globLogW') || isempty(postk.globLogW) ...
        || ~isfield(postk, 'globTable') || isempty(postk.globTable)
    [pgospa, pgospa_decomp] = PGOSPA(x_r, truth.x, truth.X, ...
        zeros(1,0), zeros(2,0), zeros(2,2,0), p, c, alpha);
    return
end

globLogW = postk.globLogW(:);
logWNorm = globLogW - logsumexp(globLogW);
weights = exp(logWNorm);

A = size(postk.globTable, 1);
pgospa = 0;
pgospa_decomp = struct('localisation', 0, 'existence_mismatch', 0, 'missed', 0, 'false', 0);

for a = 1:A
    mb = buildMBFromGlobal(postk.trackBern, postk.globTable(a, :));
    [dA, decompA] = PGOSPA(x_r, truth.x, truth.X, mb.r, mb.x, mb.X, p, c, alpha);

    wa = weights(a);
    pgospa = pgospa + wa * dA;
    pgospa_decomp.localisation = pgospa_decomp.localisation + wa * decompA.localisation;
    pgospa_decomp.existence_mismatch = pgospa_decomp.existence_mismatch + wa * decompA.existence_mismatch;
    pgospa_decomp.missed = pgospa_decomp.missed + wa * decompA.missed;
    pgospa_decomp.false = pgospa_decomp.false + wa * decompA.false;
end

end

function truth = buildTruthSet(data, k)
valid = ~isnan(squeeze(data.truthTracks(1, k, :)));
idx = find(valid);

N = numel(idx);
truth.x = zeros(2, N);
truth.X = zeros(2, 2, N);

for n = 1:N
    o = idx(n);
    truth.x(:, n) = data.truthTracks(1:2, k, o);
    truth.X(:, :, n) = data.truthExtents(:, :, k, o);
end
end

function mb = buildMBFromGlobal(trackBern, row)
K = numel(trackBern);

r = zeros(1, K);
y = zeros(2, K);
Y = zeros(2, 2, K);
count = 0;

for s = 1:K
    hIdx = row(s);
    if hIdx == 0
        continue
    end

    pool = trackBern{s};
    if hIdx > numel(pool)
        continue
    end

    hypo = pool{hIdx};
    if ~isfield(hypo, 'existence') || isempty(hypo.existence)
        continue
    end

    rs = hypo.existence;
    if rs <= 0
        continue
    end

    if ~isfield(hypo, 'ggiw') || isempty(hypo.ggiw)
        continue
    end

    g = hypo.ggiw;
    if ~isfield(g, 'mean') || ~isfield(g, 'V') || ~isfield(g, 'v')
        continue
    end
    if isempty(g.mean) || isempty(g.V) || isempty(g.v) || g.v <= 3
        continue
    end

    count = count + 1;
    r(count) = rs;
    y(:, count) = g.mean(1:2);
    Y(:, :, count) = g.V / (g.v - 3);
end

mb.r = r(1:count);
mb.x = y(:, 1:count);
mb.X = Y(:, :, 1:count);
end
