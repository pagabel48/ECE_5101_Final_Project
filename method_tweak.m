function [rxSignalTweak, est] = method_tweak(rxSignal, p)
% Our tweaked method - extends periodic update by recording the phase
% correction history and interpolating between update points so the
% compensation is smooth instead of jumping every K symbols.

    rxSignal = rxSignal(:);
    sigLen = length(rxSignal);
    symLen = p.N + p.cpLen;

    preambleSamples = 0;
    if p.useRepeatedPreamble
        preambleSamples = p.preambleLen * p.numPreambleRepeats;
    end
    payloadStart = preambleSamples + 1;

    % start with the same coarse estimate as one-shot
    rxPreamble = rxSignal(1:preambleSamples);
    epsilonCoarse = estimate_doppler_from_preamble(rxPreamble, p);

    rxCoarse = compensate_signal(rxSignal, epsilonCoarse, p);
    rxData = rxCoarse(payloadStart:end);
    rxCP = reshape(rxData, symLen, p.numSymbols);
    rxNoCP = rxCP(p.cpLen+1:end, :);

    % run the same decision-directed loop as periodic update,
    % but save the phase value every time we do an update
    K = p.updateInterval;
    phaseComp = 0;

    updateSamplePos = payloadStart;  % seed at payload start
    updatePhaseVal  = 0;

    for m = 1:p.numSymbols
        y = fft(rxNoCP(:, m), p.N);
        y_corr = y * exp(-1j * phaseComp);

        if mod(m-1, K) == 0
            x_dec = qpsk_slicer(y_corr);
            phi_hat = angle(sum(y_corr .* conj(x_dec)));
            phaseComp = phaseComp + phi_hat;

            % log sample position and accumulated phase at this update
            spos = payloadStart + (m-1)*symLen + round(symLen/2);
            updateSamplePos(end+1) = spos;
            updatePhaseVal(end+1)  = phaseComp;
        end
    end

    % interpolate the phase history to every sample -
    % this makes compensation smooth rather than step-wise
    sampleAxis = (1:sigLen).';
    phasePerSample = interp1(updateSamplePos(:), updatePhaseVal(:), ...
                             sampleAxis, 'linear', 'extrap');

    n = (0:sigLen-1).';
    totalPhase = 2*pi * epsilonCoarse * n / p.N + phasePerSample;
    rxSignalTweak = rxSignal .* exp(-1j * totalPhase);

    est.epsilonCoarse     = epsilonCoarse;
    est.epsilonRefined    = epsilonCoarse;
    est.epsilonApproxMean = epsilonCoarse;
    est.updateSamplePos   = updateSamplePos;
    est.updatePhaseVal    = updatePhaseVal;
    est.phasePerSample    = phasePerSample;
    est.epsilonGrid       = [];
    est.score             = [];
end

function x_dec = qpsk_slicer(y)
    x_dec = (sign(real(y)) + 1j*sign(imag(y))) / sqrt(2);
end

function epsilon_hat = estimate_doppler_from_preamble(rxPreamble, p)
    L = p.preambleLen;
    r1 = rxPreamble(1:L);
    r2 = rxPreamble(L+1:2*L);
    P = sum(r2 .* conj(r1));
    epsilon_hat = angle(P) * p.N / (2*pi*L);
end