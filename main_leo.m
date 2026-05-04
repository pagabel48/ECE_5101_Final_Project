clc; clear; close all;

p = params();

% LEO: faster and larger Doppler variation than GEO
p.epsilon_dc     = 0.05;
p.epsilon_amp    = 0.025;
p.epsilon_cycles = 4.0;
p.plotConstellation = true;

[txSignal, txBits, txMatrix, preamble] = ofdm_tx(p);

x = linspace(0, 1, length(txSignal)).';
epsilonVec = p.epsilon_dc + p.epsilon_amp * sin(2*pi*p.epsilon_cycles*x);
phase = cumsum(2*pi * epsilonVec / p.N);
txSignalDoppler = txSignal .* exp(1j * phase);

signalPower = mean(abs(txSignalDoppler).^2);
noisePower  = signalPower / 10^(p.snr_dB/10);
noise = sqrt(noisePower/2) * (randn(size(txSignalDoppler)) + 1j*randn(size(txSignalDoppler)));
rxSignal = txSignalDoppler + noise;

figure;
plot(epsilonVec, 'LineWidth', 1.5); grid on;
title('LEO Doppler Profile');
xlabel('Sample Index'); ylabel('Normalized Epsilon');

%% Method 0: No Compensation
[rxBits, rxSymbols, ~] = method_no_comp(rxSignal, p);
ber = mean(txBits ~= rxBits);
fprintf('BER (no compensation)  = %g\n', ber);
plot_constellation(rxSymbols, 'No Compensation (LEO)', p);

%% Method 1: One-Shot
[rxBits1, rxSymbols1, ~, ~] = method_one_shot(rxSignal, p);
ber1 = mean(txBits ~= rxBits1);
fprintf('BER (one-shot)         = %g\n', ber1);
plot_constellation(rxSymbols1, 'One-Shot Compensation (LEO)', p);

%% Method 2: Perfect
rxSignalPerfect = compensate_signal_vec(rxSignal, epsilonVec, p);
[rxBitsPerfect, rxSymbolsPerfect, ~] = method_no_comp(rxSignalPerfect, p);
berPerfect = mean(txBits ~= rxBitsPerfect);
fprintf('BER (perfect)          = %g\n', berPerfect);
plot_constellation(rxSymbolsPerfect, 'Perfect Compensation (LEO)', p);

%% Method 3: Periodic Update
[rxBits2, rxSymbols2, ~, ~, ~] = method_periodic_update(rxSignal, p);
ber2 = mean(txBits ~= rxBits2);
fprintf('BER (periodic update)  = %g\n', ber2);
plot_constellation(rxSymbols2, 'Periodic Update (LEO)', p);

%% Method 4: Tweaked
[rxSignalTweak, ~] = method_tweak(rxSignal, p);
[rxBitsTweak, rxSymbolsTweak, ~] = method_no_comp(rxSignalTweak, p);
berTweak = mean(txBits ~= rxBitsTweak);
fprintf('BER (tweaked)          = %g\n', berTweak);
plot_constellation(rxSymbolsTweak, 'Tweaked Method (LEO)', p);

fprintf('\n=========== BER Summary (LEO) ===========\n');
fprintf('No Compensation  : %g\n', ber);
fprintf('One-Shot         : %g\n', ber1);
fprintf('Perfect          : %g\n', berPerfect);
fprintf('Periodic Update  : %g\n', ber2);
fprintf('Tweaked Method   : %g\n', berTweak);
fprintf('=========================================\n');

function plot_constellation(symbols, titleStr, p)
    if ~p.plotConstellation; return; end
    figure;
    plot(symbols, '.', 'MarkerSize', 3);
    title(titleStr); xlabel('In-Phase'); ylabel('Quadrature'); grid on;
end

function rxComp = compensate_signal_vec(rxSignal, epsilonVec, p)
    phase  = cumsum(2*pi * epsilonVec(:) / p.N);
    rxComp = rxSignal(:) .* exp(-1j * phase);
end