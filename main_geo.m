clc; clear; close all;

p = params();

% GEO-like: slow sinusoidal Doppler drift
p.epsilon_dc     = 0.025;
p.epsilon_amp    = 0.006;
p.epsilon_cycles = 0.25;

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
title('GEO-like Doppler Profile');
xlabel('Sample Index'); ylabel('Normalized Epsilon');

%% No Compensation
[rxBits, rxSymbols, ~] = method_no_comp(rxSignal, p);
ber = mean(txBits ~= rxBits);
fprintf('BER (no compensation) = %g\n', ber);
if p.plotConstellation
    figure; plot(rxSymbols, '.'); grid on;
    title('No Compensation (GEO-like)'); xlabel('In-Phase'); ylabel('Quadrature');
end

%% One-Shot
[rxBits1, rxSymbols1, ~, epsilon_hat] = method_one_shot(rxSignal, p);
ber1 = mean(txBits ~= rxBits1);
fprintf('Estimated epsilon (one-shot) = %g\n', epsilon_hat);
fprintf('BER (one-shot) = %g\n', ber1);
if p.plotConstellation
    figure; plot(rxSymbols1, '.'); grid on;
    title('One-Shot Compensation (GEO-like)'); xlabel('In-Phase'); ylabel('Quadrature');
end

%% Perfect Compensation
rxSignalPerfect = compensate_signal_vec(rxSignal, epsilonVec, p);
[rxBitsPerfect, rxSymbolsPerfect, ~] = method_no_comp(rxSignalPerfect, p);
berPerfect = mean(txBits ~= rxBitsPerfect);
fprintf('BER (perfect) = %g\n', berPerfect);
if p.plotConstellation
    figure; plot(rxSymbolsPerfect, '.'); grid on;
    title('Perfect Compensation (GEO-like)'); xlabel('In-Phase'); ylabel('Quadrature');
end

%% Periodic Update
[rxBits2, rxSymbols2, ~, ~, ~] = method_periodic_update(rxSignal, p);
ber2 = mean(txBits ~= rxBits2);
fprintf('BER (periodic update) = %g\n', ber2);
if p.plotConstellation
    figure; plot(rxSymbols2, '.'); grid on;
    title('Periodic Update (GEO-like)'); xlabel('In-Phase'); ylabel('Quadrature');
end

%% Tweaked Method
[rxSignalTweak, ~] = method_tweak(rxSignal, p);
[rxBitsTweak, rxSymbolsTweak, ~] = method_no_comp(rxSignalTweak, p);
berTweak = mean(txBits ~= rxBitsTweak);
fprintf('BER (tweaked method) = %g\n', berTweak);
if p.plotConstellation
    figure; plot(rxSymbolsTweak, '.'); grid on;
    title('Tweaked Method (GEO-like)'); xlabel('In-Phase'); ylabel('Quadrature');
end

fprintf('\n=========== BER Summary (GEO-like) ===========\n');
fprintf('No Compensation  : %g\n', ber);
fprintf('One-Shot         : %g\n', ber1);
fprintf('Perfect          : %g\n', berPerfect);
fprintf('Periodic Update  : %g\n', ber2);
fprintf('Tweaked Method   : %g\n', berTweak);
fprintf('===============================================\n');

function rxComp = compensate_signal_vec(rxSignal, epsilonVec, p)
    phase  = cumsum(2*pi * epsilonVec(:) / p.N);
    rxComp = rxSignal(:) .* exp(-1j * phase);
end