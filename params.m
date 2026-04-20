function p = params()
%PARAMS  Parameters for OFDM Doppler compensation project
    %  OFDM system parameters
    p.N               = 64; 
    p.cpLen           = 16;
    p.numSymbols      = 100;
    p.M               = 4;
    p.bitsPerSymbol   = log2(p.M);

    %  Preamble
    p.preambleLen     = 64;
    p.useRepeatedPreamble = true;
    p.numPreambleRepeats  = 2;

    %  Channel
    p.snr_dB          = 20;
    p.channelType     = 'awgn'; 

    %  Doppler
    p.dopplerMode     = 'constant'; 

    % Normalized Doppler
    p.epsilon0        = 0.05;    % initial normalized frequency offset
    p.epsilonSlope    = 1e-4;    % slope for linear/ time-varying cases
    p.epsilonNoiseStd = 0.0;     % optional small random fluctuation

    %  Periodic update settings
    p.updateInterval  = 5;       % update every K OFDM symbols
    p.historyLength   = 5;       % how many past estimates to keep
    p.tweakMode       = 'moving_average';

    %  Simulation options
    p.randomSeed      = 42;
    p.plotConstellation = true;
    p.verbose         = true;

    %  Derived quantities
    p.numDataBits = p.N * p.numSymbols * p.bitsPerSymbol;
    p.symbolLen   = p.N + p.cpLen;   % time-domain samples per OFDM symbol

end