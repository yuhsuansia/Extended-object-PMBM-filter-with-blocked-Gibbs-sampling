function bestGlob = selectBestGlobalHypo(globTable, globLogW, numTracks)
%SELECTBESTGLOBALHYPO Select the highest-weight predicted global hypothesis.
if ~isempty(globLogW)
    [~, idxMax] = max(globLogW);
    bestGlob = globTable(idxMax, 1:numTracks);
else
    bestGlob = uint32(ones(1, numTracks));
end
end
