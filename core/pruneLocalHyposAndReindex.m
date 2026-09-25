function [trackBernUpd, globTable, globLogW] = pruneLocalHyposAndReindex(trackBernUpd, globTable, globLogW, existThr)
%PRUNELOCALHYPOSANDREINDEX Prune unused or weak local hypotheses and reindex globals.
K = numel(trackBernUpd);

if ~isempty(globTable)
    for s = 1:K
        pool = trackBernUpd{s};
        L = numel(pool);
        if L == 0, continue; end
 
        used = false(L, 1);
        col = globTable(:, s);
        col = col(col > 0);
        if ~isempty(col)
            used(unique(double(col))) = true;
        end

        if any(~used)
            for h = find(~used).'
                pool{h}.existence = 0;
            end
            trackBernUpd{s} = pool;
        end
    end
end

for s = 1:K
    pool = trackBernUpd{s};
    L = numel(pool);

    r = zeros(L,1);
    for h = 1:L
        r(h) = pool{h}.existence; 
    end
    keepH = (r >= existThr);
    nKeep = nnz(keepH);

    if nKeep == 0
        trackBernUpd{s} = {};
        globTable(:,s)  = 0;
        continue
    end

    old2new = zeros(L,1,'uint32');
    old2new(keepH) = uint32(1:nKeep);

    col = globTable(:,s);
    mask = (col > 0);
    if any(mask)
        col(mask) = old2new(col(mask));
        globTable(:,s) = col;
    end

    trackBernUpd{s} = pool(keepH);
end

nonEmpty = ~cellfun(@isempty, trackBernUpd);
trackBernUpd = trackBernUpd(nonEmpty);
globTable    = globTable(:, nonEmpty);

if ~isempty(globTable)
    [U, ~, ic] = unique(globTable, 'rows', 'stable');
    merged = -inf(size(U,1),1);
    for g = 1:size(U,1)
        merged(g) = logsumexp(globLogW(ic == g));
    end
    globTable = U;
    globLogW  = merged;
end

if ~isempty(globLogW)
    globLogW = globLogW - logsumexp(globLogW);
end
end
