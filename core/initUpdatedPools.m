function [trackBernUpd, hypoCache] = initUpdatedPools(numSlots)
%INITUPDATEDPOOLS Allocate local hypothesis pools and per-slot lookup maps.
trackBernUpd = cell(numSlots, 1);
hypoCache    = cell(numSlots, 1);
for s = 1:numSlots
    hypoCache{s} = containers.Map('KeyType','char','ValueType','uint32');
end
end
