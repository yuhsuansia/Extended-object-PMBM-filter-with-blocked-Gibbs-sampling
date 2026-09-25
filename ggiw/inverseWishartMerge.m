function [mergedv, mergedV] = inverseWishartMerge(vs, Vs, weights)
%INVERSEWISHARTMERGE Merge extent factors while preserving their weighted mean.
d = size(Vs,1);

vs      = double(vs(:));
weights = double(weights(:));
weights = weights / sum(weights);

% pick merged dof (simple choice)
mergedv = sum(weights .* vs);
mergedv = max(mergedv, d + 1 + 1e-6);  % ensure mean exists

% match mean: E[X] = V/(v-d-1)
M = zeros(d,d);
for i = 1:numel(vs)
    M = M + weights(i) * (Vs(:,:,i) / (vs(i) - d - 1));
end
M = 0.5*(M+M.');

mergedV = (mergedv - d - 1) * M;
mergedV = 0.5*(mergedV+mergedV.');
end
