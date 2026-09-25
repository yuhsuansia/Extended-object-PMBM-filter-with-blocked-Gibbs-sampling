function [idx, found] = getLocalHypoIdx(cache, key)
% Return existing index if present, without building the hypothesis.
if isKey(cache, key)
    idx = cache(key);
    found = true;
else
    idx = uint32(0);
    found = false;
end
end