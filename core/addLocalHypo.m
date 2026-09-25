function [pool, cache, idx] = addLocalHypo(pool, cache, key, hypo)
% Add hypothesis assuming key is not present.
pool{end+1,1} = hypo;
idx = uint32(numel(pool));
cache(key) = idx;
end