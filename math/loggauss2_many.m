function ll = loggauss2_many(z, muMat, S11, S12, S22)
%LOGGAUSS2_MANY Evaluate 2-D Gaussian log densities for one or more candidates.
dz1 = z(1) - muMat(1,:).';
dz2 = z(2) - muMat(2,:).';

detS = S11.*S22 - S12.*S12;
inv11 =  S22 ./ detS;
inv22 =  S11 ./ detS;
inv12 = -S12 ./ detS;

quad = inv11.*(dz1.^2) + 2*inv12.*(dz1.*dz2) + inv22.*(dz2.^2);

ll = -0.5*(2*log(2*pi) + log(detS) + quad);
end
