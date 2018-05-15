%================================================
%   This program reduces the low-fidelity state-
%   space matrices from FindF16Dynamics.m to the
%   a reduced form and saves them for later use.    
%   Diffence to the original (openloop analysis)
%   version is the inclusion of the altitude h 
%   as a state.
%   Also, altitude is converted to m and speed
%   to m/s
%
%   Author: Lukas
% 
%================================================
clear;


% load low-fidelity state space matrices created by FindF16Dynamics.m
load sep_lofi_ss_terrainfollow

% Establish conversion factor from feet to meters
ft2m = 0.3048;

% -------------------longitudinal----------------

% remove engine and actuator dynamics
A_lon_ac = A_longitude_lo(1:5,1:5);
B_lon_ac = A_longitude_lo(1:5,[6 7]);

% remove altitude and actuator state from output matrices
C_lon_ac = C_longitude_lo(1:5,1:5);
D_lon_ac = zeros(5,2);

% convert ft to m
A_lon_ac(:,1) = ft2m*A_lon_ac(:,1); % height in ft to m
A_lon_ac(:,3) = ft2m*A_lon_ac(:,3); % speed in ft/s to m/s

% --------------rename and save------------------
A = A_lon_ac;
B = B_lon_ac;
C = C_lon_ac;
D = D_lon_ac;

save redu_ss_terrainfollow A B C D



