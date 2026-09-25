function idx = sampleLogCategorical(logw, skipSingletonDraw)
%SAMPLELOGCATEGORICAL Draw an index from unnormalized log weights.
% Subset sampling skips the RNG draw for a single finite candidate.
if nargin < 2
    skipSingletonDraw = false;
end
logw = logw(:);
if isempty(logw)
    error('sampleLogCategorical:EmptyInput', 'At least one log weight is required.');
end
if any(isnan(logw))
    error('sampleLogCategorical:NaNWeights', 'Log weights must not contain NaN values.');
end

m = max(logw);
if m == -inf
    error('sampleLogCategorical:AllImpossible', 'At least one log weight must be finite.');
end
if m == inf
    dom = find(logw == inf);
    idx = dom(randi(numel(dom)));
    return
end

if skipSingletonDraw && isscalar(logw)
    idx = 1;
    return
end

w = exp(logw - m);
c = cumsum(w);
u = rand * c(end);
idx = find(u <= c, 1, 'first');
if isempty(idx)
    idx = numel(logw);
end
end
