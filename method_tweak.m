function [rxSignal_tweak, est] = method_tweak(rxSignal, p)
%METHOD_TWEAK
% Tweaked method for current constant-epsilon model:
% one-shot estimate + tighter local constellation refinement

    rxSignal = rxSignal(:);

    %% Part 1: coarse epsilon from one-shot
    [~, ~, ~, epsilonCoarse] = method_one_shot(rxSignal, p);

    %% Part 2: tighter local refinement
    searchSpan = 0.03;     % much tighter than before
    numGrid = 101;

    epsilonGrid = linspace(epsilonCoarse - searchSpan, ...
                           epsilonCoarse + searchSpan, ...
                           numGrid);

    score = zeros(numGrid,1);

    for k = 1:numGrid
        y = compensate_signal(rxSignal, epsilonGrid(k), p);
        [~, rxSymbolsCand, ~] = method_no_comp(y, p);
        score(k) = constellation_score_qpsk(rxSymbolsCand);
    end

    [~, bestIdx] = min(score);
    epsilonRefined = epsilonGrid(bestIdx);

    %% Part 3: final compensation
    rxSignal_tweak = compensate_signal(rxSignal, epsilonRefined, p);

    %% Save debug info
    est.epsilonCoarse = epsilonCoarse;
    est.epsilonRefined = epsilonRefined;
    est.epsilonGrid = epsilonGrid;
    est.score = score;
end


function score = constellation_score_qpsk(rxSymbols)
% Smaller score is better

    rxSymbols = rxSymbols(:);

    ref = [ 1+1j;
            1-1j;
           -1+1j;
           -1-1j ] / sqrt(2);

    dist2 = zeros(length(rxSymbols), 4);

    for m = 1:4
        dist2(:,m) = abs(rxSymbols - ref(m)).^2;
    end

    minDist2 = min(dist2, [], 2);
    score = mean(minDist2);
end