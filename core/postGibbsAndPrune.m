function [trackBernUpd, globTableUpd, globLogWUpd] = postGibbsAndPrune(trackBernUpd, globTableUpd, par)
%POSTGIBBSANDPRUNE Weight unique sampled hypotheses and prune the posterior.
globTableUpd = unique(globTableUpd, 'rows', 'stable');
Aupd = size(globTableUpd,1);

numSlots = size(globTableUpd,2);

logWslot = cell(numSlots,1);
for s = 1:numSlots
    pool = trackBernUpd{s};
    if isempty(pool)
        logWslot{s} = [];
        continue
    end
    Ls = numel(pool);
    lw = zeros(Ls,1);
    for h = 1:Ls
        lw(h) = pool{h}.logWeight;
    end
    logWslot{s} = lw;
end

globLogWUpd = zeros(Aupd,1);
for a = 1:Aupd
    row = globTableUpd(a,:);
    acc = 0;
    for s = 1:numSlots
        h = row(s);
        if h > 0
            acc = acc + logWslot{s}(h);
        end
    end
    globLogWUpd(a) = acc;
end

globLogWUpd = globLogWUpd - logsumexp(globLogWUpd);

keep = (globLogWUpd >= log(par.globalHypoPruneThreshold));
if any(keep)
    globTableUpd = globTableUpd(keep,:);
    globLogWUpd  = globLogWUpd(keep);
end

[trackBernUpd, globTableUpd, globLogWUpd] = ...
    pruneLocalHyposAndReindex(trackBernUpd, globTableUpd, globLogWUpd, par.localHypoExistenceThreshold);
end
