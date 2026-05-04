function [rxBits, rxSymbols, rxMatrix, epsilon_hat, phaseHist] = method_periodic_update(rxSignal, p)
    preambleSamples = 0;
    if p.useRepeatedPreamble
        preambleSamples = p.preambleLen * p.numPreambleRepeats;
    end

    % coarse one-shot estimate first
    rxPreamble  = rxSignal(1:preambleSamples);
    epsilon_hat = estimate_doppler_from_preamble(rxPreamble, p);

    rxSignalComp = compensate_signal(rxSignal, epsilon_hat, p);
    rxData = rxSignalComp(preambleSamples + 1:end);
    rxCP   = reshape(rxData, p.N + p.cpLen, p.numSymbols);
    rxNoCP = rxCP(p.cpLen + 1:end, :);

    rxMatrix  = zeros(p.N, p.numSymbols);
    phaseHist = zeros(p.numSymbols, 1);
    phaseComp = 0;

    for m = 1:p.numSymbols
        y = fft(rxNoCP(:, m), p.N);
        y = y * exp(-1j * phaseComp);

        % re-estimate residual phase every K symbols
        if mod(m-1, p.updateInterval) == 0
            x_dec   = qpsk_slicer(y);
            phi_hat = angle(sum(y .* conj(x_dec)));
            phaseComp = phaseComp + phi_hat;
            y = y * exp(-1j * phi_hat);
        end

        rxMatrix(:, m) = y;
        phaseHist(m)   = phaseComp;
    end

    rxSymbols = rxMatrix(:);
    rxBits = zeros(length(rxSymbols)*2, 1);
    rxBits(1:2:end) = real(rxSymbols) > 0;
    rxBits(2:2:end) = imag(rxSymbols) > 0;
end

function x_dec = qpsk_slicer(y)
    x_dec = sign(real(y)) + 1j*sign(imag(y));
    x_dec = x_dec / sqrt(2);
    x_dec(real(x_dec) == 0) = 1/sqrt(2) + 1j*imag(x_dec(real(x_dec) == 0));
    x_dec(imag(x_dec) == 0) = real(x_dec(imag(x_dec) == 0)) + 1j/sqrt(2);
end

function epsilon_hat = estimate_doppler_from_preamble(rxPreamble, p)
    L  = p.preambleLen;
    r1 = rxPreamble(1:L);
    r2 = rxPreamble(L+1:2*L);
    P  = sum(r2 .* conj(r1));
    epsilon_hat = angle(P) * p.N / (2*pi*L);
end