function [rxBits, rxSymbols, rxMatrix, epsilon_hat] = method_one_shot(rxSignal, p)
    % preamble length in samples
    preambleSamples = 0;
    if p.useRepeatedPreamble
        preambleSamples = p.preambleLen * p.numPreambleRepeats;
    end

    % one shot Doppler
    rxPreamble = rxSignal(1:preambleSamples);
    epsilon_hat = estimate_doppler_from_preamble(rxPreamble, p);

    % compensate received signal
    rxSignalComp = compensate_signal(rxSignal, epsilon_hat, p);

    % remove preamble and CP
    rxData = rxSignalComp(preambleSamples + 1:end);
    rxCP = reshape(rxData, p.N + p.cpLen, p.numSymbols);
    rxNoCP = rxCP(p.cpLen + 1:end, :);

    % FFT and QPSK hard decision
    rxMatrix = fft(rxNoCP, p.N, 1);
    rxSymbols = rxMatrix(:);
    rxBits = zeros(length(rxSymbols)*2, 1);
    rxBits(1:2:end) = real(rxSymbols) > 0;
    rxBits(2:2:end) = imag(rxSymbols) > 0;
end

function epsilon_hat = estimate_doppler_from_preamble(rxPreamble, p)
    % split repeated preamble
    L = p.preambleLen;
    r1 = rxPreamble(1:L);
    r2 = rxPreamble(L+1:2*L);
    % phase-based offset estimate
    P = sum(r2 .* conj(r1));
    phi = angle(P);

    epsilon_hat = phi * p.N / (2*pi*L);
end