function poissonUndetected = predictPPP(poissonUndetected, birth, par, cache)
%PREDICTPPP Update/prune undetected intensity, predict, and append births.
poissonUndetected = pppMissedDetection(poissonUndetected);

if ~isempty(poissonUndetected)
    keepPPP = [poissonUndetected.logWeight] > log(par.poissonUndetectedPruneThresh);
    poissonUndetected = poissonUndetected(keepPPP);
end

for i = 1:numel(poissonUndetected)
    poissonUndetected(i).logWeight = poissonUndetected(i).logWeight + cache.logPS;
    poissonUndetected(i).ggiw      = ggiwPrediction(poissonUndetected(i).ggiw, par, cache);
end

poissonUndetected(end+1) = birth;
end
