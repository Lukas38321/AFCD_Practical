%================================================
%   This program reduces the low-fidelity state-
%   space matrices from FindF16Dynamics.m to the
%   traditional (Flight Dynamics) form and saves
%   them for later use
%
%   Author: Lukas
% 
%================================================


% load low-fidelity state space matrices created by FindF16Dynamics.m
load sep_lofi_ss

% -------------------longitudinal----------------

% remove altitude from complete SS matrix
A_lon_OL = A_longitude_lo(2:7,2:7);

% remove engine and actuator dynamics
A_lon_ac = A_lon_OL(1:4,1:4);
B_lon_ac = A_lon_OL(1:4,6);

% remove altitude and actuator state from output matrices
C_lon_ac = C_longitude_lo(2:5,2:5);
D_lon_ac = zeros(4,1);

% -------------------lateral---------------------

% remove altitude from complete SS matrix
A_lat_OL = A_lateral_lo([1 4 5 6 7 8 9],[1 4 5 6 7 8 9]);

% remove engine and actuator dynamics
A_lat_ac = A_lat_OL(1:4,1:4);
B_lat_ac = A_lat_OL(1:4,[6 7]);

% remove altitude and actuator state from output matrices
C_lat_ac = C_lateral_lo([1 4 5 6],[1 4 5 6]);
D_lat_ac = zeros(4,2);


% --------------rename and save------------------
A_lon = A_lon_ac;
B_lon = B_lon_ac;
C_lon = C_lon_ac;
D_lon = D_lon_ac;

A_lat = A_lat_ac;
B_lat = B_lat_ac;
C_lat = C_lat_ac;
D_lat = D_lat_ac;

save redu_ss A_lon B_lon C_lon D_lon A_lat B_lat C_lat D_lat



