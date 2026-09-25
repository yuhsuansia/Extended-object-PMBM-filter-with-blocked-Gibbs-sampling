function updatedLocalHypo = bernoulliUpdate(predictedLocalHypo, measurements, parameters)
%BERNOULLIUPDATE Update an existing Bernoulli with an assigned cell.
r  = predictedLocalHypo.existence;
lw = predictedLocalHypo.logWeight;

[ggiw_upd, logLik] = ggiwUpdate(predictedLocalHypo.ggiw, measurements, parameters);

updatedLocalHypo = predictedLocalHypo;
updatedLocalHypo.ggiw = ggiw_upd;
updatedLocalHypo.existence = 1;
updatedLocalHypo.logWeight = lw + log(r) + logLik;

end
