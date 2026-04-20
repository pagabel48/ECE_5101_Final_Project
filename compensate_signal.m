function rxSignalComp = compensate_signal(rxSignal, epsilon_hat, p)
    n = (0:length(rxSignal)-1).';
    rxSignalComp = rxSignal .* exp(-1j * 2*pi * epsilon_hat * n / p.N);
end