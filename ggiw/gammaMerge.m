function [mergedAlpha, mergedBeta] = gammaMerge(alphas, betas, weights)
%GAMMAMERGE Merge gamma factors, preserving their weighted rate mean.
alphas  = alphas(:);
betas   = betas(:);
weights = weights(:);

alpha = sum(weights .* alphas);

mergedAlpha = alpha;
mergedBeta  = mergedAlpha / sum(weights .* (alphas ./ betas));
end
