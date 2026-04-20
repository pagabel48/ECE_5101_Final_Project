clc; clear; close all;


p = params();

snrVec = 5:5:30;
numSNR = length(snrVec);

oneShot = zeros(numSNR, 100);

for s = 1:numSNR

    snr = snrVec(s);

    for i = 1:100

        % -----------------------------
        % TX
        % -----------------------------
        [txSignal, txBits, txMatrix, preamble] = ofdm_tx(p);

        Nsig = length(txSignal);

        % -----------------------------
        % KA-BAND CHANNEL (DOPPLER)
        % -----------------------------
        h = generatechannel(Nsig, p);

        txSignalChan = txSignal .* h;

        % -----------------------------
        % AWGN
        % -----------------------------
        sigPower = mean(abs(txSignalChan).^2);
        snrLin = 10^(snr/10);
        noisePower = sigPower / snrLin;

        noise = sqrt(noisePower/2) * ...
            (randn(size(txSignalChan)) + 1j*randn(size(txSignalChan)));

        rxSignal = txSignalChan + noise;

        % -----------------------------
        % RX
        % -----------------------------
        [rxBits, ~, ~, ~] = method_one_shot(rxSignal, p);

        % BER
        numErr = sum(txBits ~= rxBits);
        oneShot(s, i) = numErr / length(txBits);

    end
end

% -----------------------------
% Plot BER vs SNR
% -----------------------------
plot(snrVec, mean(oneShot,2), 'o-');
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('One-shot OFDM under Ka-band Doppler');