function p = params()
% system parameters

% OFDM
p.N            = 64;
p.cpLen        = 16;
p.numSymbols   = 100;
p.M            = 4;
p.bitsPerSymbol = log2(p.M);

% preamble
p.preambleLen        = 64;
p.useRepeatedPreamble = true;
p.numPreambleRepeats  = 2;

% channel / noise
p.snr_dB     = 20;
p.channelType = 'awgn';

% Doppler
p.epsilon0        = 0.05;
p.epsilonSlope    = 1e-4;
p.epsilonNoiseStd = 0.0;

% periodic update interval
p.updateInterval = 5;
p.historyLength  = 5;
p.tweakMode      = 'moving_average';

p.randomSeed        = 42;
p.plotConstellation = true;
p.verbose           = true;

% derived
p.numDataBits = p.N * p.numSymbols * p.bitsPerSymbol;
p.symbolLen   = p.N + p.cpLen;

% GEO-like Doppler profile (used in main_geo.m)
p.channelMode    = 'geo_like_slow';
p.epsilon_dc     = 0.03;
p.epsilon_amp    = 0.008;
p.epsilon_cycles = 0.3;

end