function gwd = GaussianWassersteinDistance(m1,P1,m2,P2)
%GAUSSIANWASSERSTEINDISTANCE 2-Wasserstein distance between Gaussians.
%   gwd^2 = ||m1-m2||^2 + tr(P1 + P2 - 2*(P1^(1/2) P2 P1^(1/2))^(1/2))
%
% Assumes P1,P2 are symmetric positive (semi-)definite. Uses eigen-decomp
% instead of sqrtm for speed and stability.

% --- mean term ---
dm  = m1 - m2;
gwd = dm.'*dm;

% --- symmetrize (important for numeric stability) ---
P1 = 0.5*(P1 + P1.');
P2 = 0.5*(P2 + P2.');

% --- trace terms ---
trP1 = trace(P1);
trP2 = trace(P2);

% --- compute tr( (P1^(1/2) P2 P1^(1/2))^(1/2) ) efficiently ---
% Eigendecompose P1 = U*D*U'
[U,D] = eig(P1);
d = real(diag(D));
d(d < 0) = 0;                 % clip tiny negatives
s = sqrt(d);

% Form A = P1^(1/2) * P2 * P1^(1/2) in the eigenbasis of P1:
% P1^(1/2) = U*diag(s)*U'
% A = U * (diag(s) * (U'*P2*U) * diag(s)) * U'
B = U.' * P2 * U;
B = 0.5*(B + B.');            % symmetrize
A = (s .* B) .* s.';          % diag(s)*B*diag(s) without forming diag matrices
A = 0.5*(A + A.');

% Eigenvalues of A are the same as of U*A*U', so:
lam = real(eig(A));
lam(lam < 0) = 0;             % clip small negatives
trSqrt = sum(sqrt(lam));

% --- final GWD ---
gwd = gwd + (trP1 + trP2 - 2*trSqrt);
gwd = sqrt(gwd);
end
