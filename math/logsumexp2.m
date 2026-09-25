function s = logsumexp2(a, b)
%LOGSUMEXP2 Compute a stable log-sum-exp of two scalar values.
m = max(a, b);
if ~isfinite(m)
    s = m;
    return
end
s = m + log(exp(a - m) + exp(b - m));
end
