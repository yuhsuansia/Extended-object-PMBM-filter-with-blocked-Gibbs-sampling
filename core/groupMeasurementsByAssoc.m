function groups = groupMeasurementsByAssoc(measAssoc, numSlots)
%GROUPMEASUREMENTSBYASSOC Group measurement indices by associated slot.
% Output groups{s} contains measurement indices assigned to slot s.

groups = repmat({[]}, numSlots, 1);
numMeas = numel(measAssoc);

if numMeas == 0
    return
end

a = double(measAssoc(:));  % values in 1..numSlots
j = (1:numMeas).';

[as, ord] = sort(a);
js = j(ord);

cut = [1; find(diff(as) ~= 0) + 1; numMeas + 1];
for t = 1:numel(cut)-1
    s = as(cut(t));
    groups{s} = sort(js(cut(t):cut(t+1)-1), 'ascend').';
end
end
