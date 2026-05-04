clc; clear; close all;

snrRange  = 0:2:30;
numTrials = 30;           % average over multiple runs to smooth the curves
numSNR    = length(snrRange);

ber_geo = zeros(numSNR, 5);
ber_leo = zeros(numSNR, 5);

fprintf('Running %d SNR points x %d trials x 2 scenarios...\n', numSNR, numTrials);

for si = 1:numSNR
    snr = snrRange(si);
    fprintf('SNR = %2d dB ... ', snr);

    acc_geo = zeros(1,5);
    acc_leo = zeros(1,5);

    for trial = 1:numTrials
        % GEO
        p = params();
        p.snr_dB            = snr;
        p.plotConstellation = false;
        p.verbose           = false;
        p.randomSeed        = trial * 100 + si;
        p.epsilon_dc        = 0.03;
        p.epsilon_amp       = 0.008;
        p.epsilon_cycles    = 0.3;

        [txSig, txBits, ~, ~] = ofdm_tx(p);
        xg = linspace(0,1,length(txSig)).';
        ev = p.epsilon_dc + p.epsilon_amp*sin(2*pi*p.epsilon_cycles*xg);
        rx_geo = add_awgn(txSig .* exp(1j*cumsum(2*pi*ev/p.N)), snr);
        acc_geo = acc_geo + compute_all_bers(rx_geo, txBits, ev, p);

        % LEO
        p_leo              = p;
        p_leo.epsilon_dc   = 0.05;
        p_leo.epsilon_amp  = 0.025;
        p_leo.epsilon_cycles = 4.0;

        [txSig_l, txBits_l, ~, ~] = ofdm_tx(p_leo);
        xl = linspace(0,1,length(txSig_l)).';
        el = p_leo.epsilon_dc + p_leo.epsilon_amp*sin(2*pi*p_leo.epsilon_cycles*xl);
        rx_leo = add_awgn(txSig_l .* exp(1j*cumsum(2*pi*el/p_leo.N)), snr);
        acc_leo = acc_leo + compute_all_bers(rx_leo, txBits_l, el, p_leo);
    end

    ber_geo(si,:) = acc_geo / numTrials;
    ber_leo(si,:) = acc_leo / numTrials;
    fprintf('GEO=[%.3f %.3f %.4f %.4f %.4f]  LEO=[%.3f %.3f %.4f %.4f %.4f]\n', ...
        ber_geo(si,:), ber_leo(si,:));
end

colors  = {'#888888','#E87D2B','#3A7DC9','#E84040','#2BAF6A'};
markers = {'o','s','^','d','p'};
labels  = {'No Comp','One-Shot','Periodic Update','Tweaked','Perfect'};

figure('Position',[100 100 1100 480]);

subplot(1,2,1);
for m = 1:5
    semilogy(snrRange, max(ber_geo(:,m), 1e-5), ['-' markers{m}], ...
        'Color', colors{m}, 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', labels{m});
    hold on;
end
grid on; ylim([1e-5 1]);
xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR — GEO (slow drift)');
legend('Location','southwest','FontSize',9);

subplot(1,2,2);
for m = 1:5
    semilogy(snrRange, max(ber_leo(:,m), 1e-5), ['-' markers{m}], ...
        'Color', colors{m}, 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', labels{m});
    hold on;
end
grid on; ylim([1e-5 1]);
xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR — LEO (fast drift)');
legend('Location','southwest','FontSize',9);

sgtitle('OFDM Doppler Compensation — BER vs SNR');

%% helpers
function bers = compute_all_bers(rxSignal, txBits, epsilonVec, p)
    bers = zeros(1,5);

    [rb, ~, ~] = method_no_comp(rxSignal, p);
    bers(1) = mean(txBits ~= rb);

    [rb1, ~, ~, ~] = method_one_shot(rxSignal, p);
    bers(2) = mean(txBits ~= rb1);

    try
        [rb2, ~, ~, ~, ~] = method_periodic_update(rxSignal, p);
        bers(3) = mean(txBits ~= rb2);
    catch
        bers(3) = 0.5;
    end

    try
        [rxTweak, ~] = method_tweak(rxSignal, p);
        [rbt, ~, ~] = method_no_comp(rxTweak, p);
        bers(4) = mean(txBits ~= rbt);
    catch
        bers(4) = 0.5;
    end

    rxPerf = compensate_signal_vec(rxSignal, epsilonVec, p);
    [rbp, ~, ~] = method_no_comp(rxPerf, p);
    bers(5) = mean(txBits ~= rbp);
end

function rx = add_awgn(tx, snr_dB)
    noisePow = mean(abs(tx).^2) / 10^(snr_dB/10);
    rx = tx + sqrt(noisePow/2)*(randn(size(tx)) + 1j*randn(size(tx)));
end

function rxComp = compensate_signal_vec(rxSignal, epsilonVec, p)
    phase  = cumsum(2*pi * epsilonVec(:) / p.N);
    rxComp = rxSignal(:) .* exp(-1j * phase);
end