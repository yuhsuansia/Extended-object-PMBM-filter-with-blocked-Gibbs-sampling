function s = logsumexp(x)
%LOGSUMEXP Compute a stable log-sum-exp over all elements.
x = x(:);
if isempty(x)
    s = -inf;
    return
end
xmax = max(x);
if ~isfinite(xmax)
    s = xmax;
    return
end
s = xmax + log(sum(exp(x - xmax)));
end
