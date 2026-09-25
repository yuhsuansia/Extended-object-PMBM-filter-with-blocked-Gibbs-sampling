function cache = buildCache(par)
%BUILDCACHE Cache motion matrices, survival, and clutter constants.
cache = struct();
cache.logClutter = log(par.clutterIntensity);
cache.logPS      = log(par.pSurvive);
cache.pS         = par.pSurvive;
cache.H          = par.H;

mm = getConstantVelocityModel(par.dt, par.sigmaAcc);
cache.A = mm.A;
cache.Q = mm.Q;
end
