function updatedLocalHypo = bernoulliMissedDetection(predictedLocalHypo)
%BERNOULLIMISSEDDETECTION Update a Bernoulli with no assigned measurements.
updatedLocalHypo = predictedLocalHypo;

[updatedLocalHypo.ggiw, logLik] = ggiwMissedDetection(predictedLocalHypo.ggiw);

r = predictedLocalHypo.existence;
log1mr = log1p(-r);
t      = logLik + log(r) - log1mr;
logMD  = log1mr + log1p(exp(t));

updatedLocalHypo.logWeight = predictedLocalHypo.logWeight + logMD;
updatedLocalHypo.existence = exp(log(r) + logLik - logMD);
end
