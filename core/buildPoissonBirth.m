function birth = buildPoissonBirth(par)
%BUILDPOISSONBIRTH Construct the single-component GGIW Poisson birth intensity.
posStd = 150;
velStd = 15;
Ppos   = (posStd^2) * eye(2);
Pvel   = (velStd^2) * eye(2);
gammaScale = par.birthGammaScale;
iwDoF = par.birthIwDoF;

birth = struct();
birth.logWeight = log(par.meanBirths);
birth.ggiw = struct( ...
    'mean',  zeros(4,1), ...
    'cov',   blkdiag(Ppos, Pvel), ...
    'alpha', par.meanMeasPerObject * gammaScale, ...
    'beta',  gammaScale, ...
    'v',     iwDoF, ...
    'V',     par.extentNominal * (iwDoF - 3) );

birth.ggiw = ggiwFinalizeCache(birth.ggiw);
end
