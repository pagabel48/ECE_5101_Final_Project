clc; clear; close all;

P = params();

% returns doppler shift vector based on realistic orbit/ground station
% N  = number of samples
% P = params
fd = generatechannel(P.numDataBits, P);
plot(real(fd)); grid on;
title("Ka-band Doppler + Channel Phase Response");

function h = generatechannel(N, P)

    c  = 3e8;
    Re = 6371e3;
    mu = 3.986e14;

    fc = P.fc;
    dt = P.dt;

    lambda = c / fc;

    % LEO orbit
    altitude = 400e3 + 600e3*rand();
    r = Re + altitude;

    v = sqrt(mu / r);

    % Time
    t = (0:N-1) * dt;

    % Random geometry
    sat_init = 2*pi*rand();

    inclination = deg2rad(30 + 60*rand());

    % ground station (simplified Earth-fixed)
    gs_lat = deg2rad(-60 + 120*rand());
    gs_lon = 2*pi*rand();

    gs_vec = [cos(gs_lat)*cos(gs_lon);
          cos(gs_lat)*sin(gs_lon);
          sin(gs_lat)];

    % Output complex channel
    h = zeros(1,N);

    for k = 1:N

        theta = sat_init + (v/r)*t(k);

        % satellite position (simplified circular orbit)
        sat_vec = [cos(theta)*cos(inclination);
               sin(theta);
               sin(theta)*sin(inclination)];

        % line-of-sight
        los = sat_vec - gs_vec;
        rhat = los / norm(los);

        % velocity vector
        v_vec = v * [-sin(theta);
                  cos(theta);
                  0];

        % radial velocity
        vr = dot(v_vec, rhat);

        % Doppler shift
        fD = vr / lambda;

        % phase accumulation (THIS is key for OFDM)
        h(k) = exp(1j*2*pi*sum(fD)*dt);

    end

    % Ka-band free space loss
        d = norm(los);
    FSPL = (4*pi*d/lambda)^2;

    h = h / sqrt(FSPL);

end
