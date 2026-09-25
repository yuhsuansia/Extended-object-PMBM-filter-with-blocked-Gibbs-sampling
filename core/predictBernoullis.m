function trackBern = predictBernoullis(trackBern, par, cache)
%PREDICTBERNOULLIS Apply survival and GGIW prediction to every local hypothesis.
pS = cache.pS;
for i = 1:numel(trackBern)
    Ti = trackBern{i};
    for h = 1:numel(Ti)
        hyp = Ti{h};
        hyp.existence = pS * hyp.existence;
        hyp.ggiw      = ggiwPrediction(hyp.ggiw, par, cache);
        Ti{h} = hyp;
    end
    trackBern{i} = Ti;
end
end
