function [txSignal, txBits, txMatrix, preamble] = ofdm_tx(p)
    rng(p.randomSeed);

    txBits = randi([0 1], p.numDataBits, 1);
    bitPairs = reshape(txBits, 2, []).';
    txSymbols = ((2*bitPairs(:,1)-1) + 1j*(2*bitPairs(:,2)-1)) / sqrt(2);

    txMatrix = reshape(txSymbols, p.N, p.numSymbols);

    txIFFT = ifft(txMatrix, p.N, 1);
    txCP = [txIFFT(end-p.cpLen+1:end, :); txIFFT];

    preamble = [];
    if p.useRepeatedPreamble
        preambleFreq = exp(1j * 2*pi * rand(p.preambleLen,1));
        preambleTime = ifft(preambleFreq, p.preambleLen);
        preamble = repmat(preambleTime, p.numPreambleRepeats, 1);
    end

    txSignal = [preamble; txCP(:)];
end