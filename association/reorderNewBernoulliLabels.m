function measAssoc = reorderNewBernoulliLabels(measAssoc, numTracks, numMeas, numSlots)
%REORDERNEWBERNOULLILABELS Canonical relabeling of new-Bernoulli components.
%
% Enforces canonical labeling for new components (slots numTracks+1 ... numTracks+numMeas):
%   - If a new component is non-empty, its label becomes k = max(meas indices in that component).
%   - Hence: if any measurements assigned to slot numTracks+k, then measurement k is included,
%            and no measurement j is assigned to a new slot index < j.
%
% This does NOT modify detStates; only measAssoc is relabeled.
%
% Inputs:
%   measAssoc : (numMeas x 1) uint32/double, each entry in 1..numSlots
%   numTracks : # legacy tracks (slots 1..numTracks)
%   numMeas   : # measurements
%   numSlots  : total # slots (usually numTracks+numMeas). If omitted, uses numTracks+numMeas.

    if nargin < 4 || isempty(numSlots)
        numSlots = numTracks + numMeas;
    end

    assoc = double(measAssoc(:));

    % Mask: measurements assigned to "new slots"
    isNew = (assoc > numTracks) & (assoc <= (numTracks + numMeas));
    if ~any(isNew)
        measAssoc = uint32(assoc);
        return
    end

    % For each slot s, compute max measurement index assigned to s (0 if empty)
    maxIdxPerSlot = accumarray(assoc, (1:numMeas).', [numSlots, 1], @max, 0);

    % Build mapping oldSlot -> newSlot for new slots that are non-empty
    old2new = uint32((1:numSlots).');  % identity by default

    usedNewSlots = find( (1:numSlots).' > numTracks & maxIdxPerSlot > 0 );
    if ~isempty(usedNewSlots)
        target = uint32(numTracks + maxIdxPerSlot(usedNewSlots));

        old2new(uint32(usedNewSlots)) = target;
    end

    % Apply mapping only to measurements currently assigned to new slots
    assoc(isNew) = double(old2new(uint32(assoc(isNew))));

    measAssoc = uint32(assoc);
end