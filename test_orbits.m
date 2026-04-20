clc; clear; close all;

p = params();

subcarriers []

% Test Ka band
[txSignal, txBits, txMatrix, preamble] = ofdm_tx(p);

oneShot = zeroes(6, )

% vary SNR from 5 to 30 db

for snr = 5 : 5 : 30

% take 100 samples of procedural satellites
for i = 1:100
    [txSignal, txBits, txMatrix, preamble] = ofdm_tx(p);

    n = (0:length(txSignal)-1).';

    % time-varying Doppler (Hz per sample)
    fd = p.fd;   % must be a vector same length as txSignal

    Ts = 1 / p.fs;   % sampling period (IMPORTANT: use real fs)

    % accumulate phase from instantaneous Doppler
    phase = 2*pi * cumsum(fd) * Ts;

    txSignalDoppler = txSignal .* exp(1j * phase);

    signalPower = mean(abs(txSignalDoppler).^2);
    snrLinear = 10^(snr/10);
    noisePower = signalPower / snrLinear;

    noise = sqrt(noisePower/2) * ...
        (randn(size(txSignalDoppler)) + 1j*randn(size(txSignalDoppler)));

    rxSignal = txSignalDoppler + noise;

end
end
