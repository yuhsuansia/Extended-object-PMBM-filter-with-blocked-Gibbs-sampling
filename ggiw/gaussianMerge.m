function [mergedMean, mergedCov] = gaussianMerge(means, covs, weights)
%GAUSSIANMERGE Match the first two moments of a Gaussian mixture.
weights = weights(:);
mergedMean = means * weights;

w3 = reshape(weights, 1, 1, []);
covTerm = sum(covs .* w3, 3);

mw = means .* weights.';
secondMoment = mw * means.';

mergedCov = covTerm + secondMoment - (mergedMean * mergedMean.');
mergedCov = 0.5 * (mergedCov + mergedCov.');
end
