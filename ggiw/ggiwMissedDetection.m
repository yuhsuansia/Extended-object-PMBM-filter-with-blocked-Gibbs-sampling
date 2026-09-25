function [updatedGGIW, predictedLogLik] = ggiwMissedDetection(predictedGGIW)
%GGIWMISSEDDETECTION Condition the gamma measurement rate on an empty cell.
a = predictedGGIW.alpha;
b = predictedGGIW.beta;

predictedLogLik = -a .* log1p(1 ./ b);

updatedGGIW = predictedGGIW;
updatedGGIW.beta = b + 1;
end
