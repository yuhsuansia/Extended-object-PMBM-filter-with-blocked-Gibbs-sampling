function X = iwishrnd2_fast_fromCache(VinvL, nu)
%IWISHRND2_FAST_FROMCACHE Draw a 2-D inverse-Wishart matrix using cached factors.
t11 = sqrt(2 * randg(nu/2));
t21 = randn;
t22 = sqrt(2 * randg((nu-1)/2));

A  = VinvL * [t11 0; t21 t22];
W  = A * A.';
W  = 0.5*(W+W.');

w11 = W(1,1); w12 = W(1,2); w22 = W(2,2);
detW = w11*w22 - w12*w12;
detW = max(detW, realmin);

X = (1/detW) * [ w22  -w12;
    -w12   w11 ];
X = 0.5*(X+X.');
end
