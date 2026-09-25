function idx = sampleLogCategoricalSubset(ids, logw)
%SAMPLELOGCATEGORICALSUBSET Sample from categorical log-weights on a subset.
% ids: candidate labels; logw: corresponding log-weights.

ids = ids(:);
logw = logw(:);

if isempty(ids)
    error('sampleLogCategoricalSubset:EmptyInput', 'At least one candidate id is required.');
end
if numel(ids) ~= numel(logw)
    error('sampleLogCategoricalSubset:SizeMismatch', 'ids and logw must have the same number of elements.');
end
idx = ids(sampleLogCategorical(logw, true));
end
