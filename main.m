clc; clear; close all;

p = params();

[txSignal, txBits, txMatrix, preamble] = ofdm_tx(p);

% Apply Doppler / CFO
n = (0:length(txSignal)-1).';
epsilon = p.epsilon0;
txSignalDoppler = txSignal .* exp(1j * 2*pi * epsilon * n / p.N);

% Add AWGN
signalPower = mean(abs(txSignalDoppler).^2);
snrLinear = 10^(p.snr_dB/10);
noisePower = signalPower / snrLinear;

noise = sqrt(noisePower/2) * ...
    (randn(size(txSignalDoppler)) + 1j*randn(size(txSignalDoppler)));

rxSignal = txSignalDoppler + noise;

%% ============================================================
% Method 0: No Compensation
%% ============================================================
[rxBits, rxSymbols, rxMatrix] = method_no_comp(rxSignal, p);

numErr = sum(txBits ~= rxBits);
ber = numErr / length(txBits);

fprintf('BER (no compensation) = %g\n', ber);

if p.plotConstellation
    figure;
    plot(rxSymbols, '.');
    title('Received Constellation (No Compensation)');
    xlabel('In-Phase');
    ylabel('Quadrature');
    grid on;
end

%% ============================================================
% Method 1: One-Shot Compensation
%% ============================================================
[rxBits1, rxSymbols1, rxMatrix1, epsilon_hat] = method_one_shot(rxSignal, p);

numErr1 = sum(txBits ~= rxBits1);
ber1 = numErr1 / length(txBits);

fprintf('Estimated epsilon (one-shot) = %g\n', epsilon_hat);
fprintf('BER (one-shot compensation) = %g\n', ber1);

if p.plotConstellation
    figure;
    plot(rxSymbols1, '.');
    title('Received Constellation (One-Shot Compensation)');
    xlabel('In-Phase');
    ylabel('Quadrature');
    grid on;
end

%% ============================================================
% Method 2: Perfect Compensation
%% ============================================================
epsilon_true = p.epsilon0;
rxSignalPerfect = compensate_signal(rxSignal, epsilon_true, p);

[rxBitsPerfect, rxSymbolsPerfect, rxMatrixPerfect] = method_no_comp(rxSignalPerfect, p);

numErrPerfect = sum(txBits ~= rxBitsPerfect);
berPerfect = numErrPerfect / length(txBits);

fprintf('BER (perfect compensation) = %g\n', berPerfect);

if p.plotConstellation
    figure;
    plot(rxSymbolsPerfect, '.');
    title('Received Constellation (Perfect Compensation)');
    xlabel('In-Phase');
    ylabel('Quadrature');
    grid on;
end

%% ============================================================
% Method 3: Periodic Update
%% ============================================================
[rxBits2, rxSymbols2, rxMatrix2, epsilon_hat2, phaseHist] = method_periodic_update(rxSignal, p);

numErr2 = sum(txBits ~= rxBits2);
ber2 = numErr2 / length(txBits);

fprintf('BER (periodic update) = %g\n', ber2);

if p.plotConstellation
    figure;
    plot(rxSymbols2, '.');
    title('Received Constellation (Periodic Update)');
    xlabel('In-Phase');
    ylabel('Quadrature');
    grid on;
end

%% ============================================================
% Method 4: Tweaked Method
% Coarse estimate from repeated preamble + local refinement search
%% ============================================================
[rxSignalTweak, estTweak] = method_tweak(rxSignal, p);

% After compensation, reuse the normal receiver
[rxBitsTweak, rxSymbolsTweak, rxMatrixTweak] = method_no_comp(rxSignalTweak, p);

numErrTweak = sum(txBits ~= rxBitsTweak);
berTweak = numErrTweak / length(txBits);

epsilon_hat_tweak_coarse = estTweak.epsilonCoarse;
epsilon_hat_tweak_total  = estTweak.epsilonRefined;

fprintf('Estimated epsilon (tweak coarse) = %g\n', epsilon_hat_tweak_coarse);
fprintf('Estimated epsilon (tweak refined) = %g\n', epsilon_hat_tweak_total);
fprintf('BER (tweaked method) = %g\n', berTweak);

if p.plotConstellation
    figure;
    plot(rxSymbolsTweak, '.');
    title('Received Constellation (Tweaked Method)');
    xlabel('In-Phase');
    ylabel('Quadrature');
    grid on;
end

%% ============================================================
% Optional debug plots for tweaked method
%% ============================================================
if isfield(p, 'plotTweakDebug') && p.plotTweakDebug
    figure;
    plot(estTweak.epsilonGrid, estTweak.score, 'LineWidth', 1.5);
    grid on;
    title('Tweaked Method: Refinement Search Metric');
    xlabel('Candidate Epsilon');
    ylabel('Metric');

    hold on;
    xline(estTweak.epsilonCoarse, '--r', 'Coarse');
    xline(estTweak.epsilonRefined, '--g', 'Refined');
    hold off;
end

%% ============================================================
% Summary
%% ============================================================
fprintf('\n=========== BER Summary ===========\n');
fprintf('No Compensation       : %g\n', ber);
fprintf('One-Shot Compensation : %g\n', ber1);
fprintf('Perfect Compensation  : %g\n', berPerfect);
fprintf('Periodic Update       : %g\n', ber2);
fprintf('Tweaked Method        : %g\n', berTweak);
fprintf('===================================\n');