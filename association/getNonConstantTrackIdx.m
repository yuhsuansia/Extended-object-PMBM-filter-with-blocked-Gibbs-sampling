function actIdx = getNonConstantTrackIdx(globTable, numTracks)
%GETNONCONSTANTTRACKIDX Track indices whose global-hypothesis column varies.

activeTracks = true(1, numTracks);
for i = 1:numTracks
    col = globTable(:, i);
    if min(col) == max(col)
        activeTracks(i) = false;
    end
end

actIdx = find(activeTracks);
end
