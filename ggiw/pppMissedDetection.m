function updatedPPP = pppMissedDetection(predictedPPP)
%PPPMISSEDDETECTION Update each undetected PPP component for an empty cell.
numGGIW = numel(predictedPPP);
updatedPPP = predictedPPP;

for i = 1:numGGIW
    comp = updatedPPP(i);
    [comp.ggiw, logLik] = ggiwMissedDetection(comp.ggiw);
    comp.logWeight = comp.logWeight + logLik;
    updatedPPP(i) = comp;
end
end
