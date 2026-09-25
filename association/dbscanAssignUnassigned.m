function measAssoc = dbscanAssignUnassigned(Z, measAssoc, numTracks, par)
%DBSCANASSIGNUNASSIGNED Cluster free measurements and assign canonical new-object slots.
unassigned = find(measAssoc == 0);
if isempty(unassigned), return; end

clusterIdx = dbscan(Z(:,unassigned)', par.dbscanEps, par.dbscanMinPts);

isNoise = (clusterIdx == -1);
measAssoc(unassigned(isNoise)) = uint32(numTracks + unassigned(isNoise));

isValid = (clusterIdx > 0);
if any(isValid)
    cl  = clusterIdx(isValid);
    mid = unassigned(isValid);
    maxIdPerCluster = accumarray(cl, mid, [], @max);
    measAssoc(mid)  = uint32(numTracks + maxIdPerCluster(cl));
end
end
