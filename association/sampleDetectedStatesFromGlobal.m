function states = sampleDetectedStatesFromGlobal(trackBernUpd, globalHypo, collapsed)
%SAMPLEDETECTEDSTATESFROMGLOBAL Draw GGIW states for the current hypothesis.
% Standard Gibbs samples existence; collapsed Gibbs retains its probability.
numSlots = numel(trackBernUpd);
states = cell(numSlots, 1);
emptyState = struct('existence', 0, 'kinematicState', [], ...
    'poissonRate', [], 'extentState', []);

for s = 1:numSlots
    hIdx = globalHypo(s);
    if collapsed
        obj = emptyState;
        if hIdx == 0
            states{s} = obj;
            continue
        end
        lh = trackBernUpd{s}{hIdx};
        obj.existence = lh.existence;
        drawState = obj.existence > 0;
    else
        if hIdx == 0, continue; end
        lh = trackBernUpd{s}{hIdx};
        drawState = lh.existence > rand;
        if ~drawState, continue; end
        obj = struct();
    end

    if drawState
        g = lh.ggiw;
        obj.kinematicState = g.mean + g.covL * randn(numel(g.mean), 1);
        obj.poissonRate = randg(g.alpha) / g.beta;
        obj.extentState = iwishrnd2_fast_fromCache(g.VinvL, g.v);
    end
    states{s} = obj;
end
end
