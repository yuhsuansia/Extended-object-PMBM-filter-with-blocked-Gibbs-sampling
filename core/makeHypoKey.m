function key = makeHypoKey(parentalIdx, measIdx, isSorted)
%MAKEHYPOKEY Build a canonical string key for containers.Map (KeyType='char').

if nargin < 3
    isSorted = false;
end

% Ensure canonical ordering when needed.
if ~isempty(measIdx) && ~isSorted
    measIdx = sort(measIdx);
end

p = uint32(parentalIdx);
n = uint32(numel(measIdx));

if n == 0
    % Small, unambiguous key for missed-detection / empty group
    key = sprintf('p%u|n0', p);
    return
end

% Build "p<parent>|n<count>|m<id1>,<id2>,...,<idn>" in one pass.
fmt = ['p%u|n%u|m' repmat('%u,', 1, double(n))];
key = sprintf(fmt, p, n, uint32(measIdx));
key(end) = [];
end
