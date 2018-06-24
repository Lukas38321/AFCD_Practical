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

% ----------select parts of full matrix----------------

% remove engine and actuator dynamics
A = A_longitude_lo(1:5,1:5);
B = A_longitude_lo(1:5,[6 7]);

% remove actuator states from output matrices
C = eye(5);
D = zeros(5,2);

% -------------------save to file----------------------

save redu_ss_terrainfollow A B C D



