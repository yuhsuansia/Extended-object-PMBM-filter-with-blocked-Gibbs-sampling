function A = ggiwExtentTransform(Xhat, S)
%GGIWEXTENTTRANSFORM Principal-square-root map for the 2D extent innovation.
% N = (A*eps)*(A*eps)' with A = Xhat^(1/2)*S^(-1/2).
% Principal symmetric roots are intentional: arbitrary Cholesky factors do
% not give the same innovation for noncommuting covariance matrices.
A = rootSPD2(Xhat) / rootSPD2(S);
end

function R = rootSPD2(M)
% Closed-form principal root of a 2-by-2 SPD matrix avoids general sqrtm.
M = 0.5*(M + M.');
a = M(1,1);
b = M(1,2);
c = M(2,2);
determinant = a*c-b*b;
if a > 0 && c > 0 && determinant > 0 && isfinite(determinant)
    s = sqrt(determinant);
    R = (M + s*eye(2)) / sqrt(a+c+2*s);
else
    % Fall back when forming the determinant loses precision or overflows.
    [U, D] = eig(M);
    values = diag(D);
    assert(all(isfinite(values) & values > 0), ...
        'ggiwExtentTransform:NotSPD', 'Expected an SPD covariance.');
    R = U*diag(sqrt(values))*U.';
    R = 0.5*(R + R.');
end
end
