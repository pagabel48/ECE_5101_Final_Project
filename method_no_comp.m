function [rxBits, rxSymbols, rxMatrix] = method_no_comp(rxSignal, p)
    % preamble length in samples
    preambleSamples = 0;
    if p.useRepeatedPreamble
        preambleSamples = p.preambleLen * p.numPreambleRepeats;
    end

    % remove preamble
    rxData = rxSignal(preambleSamples + 1:end);

    % remove CP
    rxCP = reshape(rxData, p.N + p.cpLen, p.numSymbols);
    rxNoCP = rxCP(p.cpLen + 1:end, :);

    % FFT and QPSK hard decision
    rxMatrix = fft(rxNoCP, p.N, 1);
    rxSymbols = rxMatrix(:);
    rxBits = zeros(length(rxSymbols)*2, 1);
    rxBits(1:2:end) = real(rxSymbols) > 0;
    rxBits(2:2:end) = imag(rxSymbols) > 0;
end