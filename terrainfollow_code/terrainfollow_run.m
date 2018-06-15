%==========================================================================
%   This program initialises the matrices and constants required for the
%   F16 terrain follow system model and runs the Simulink model
% 
%   Author: Lukas
% 
%==========================================================================
clear;


% load reduced longitudinal state-space matrices and trim state

load redu_ss_terrainfollow
load trim_lo

% Unit conversion factor
ft2m      = 0.3048;

% simulation parameters

h0  = 1500;  % initial altitude [m]
v0  = ft2m*trim_state_lo(7); % initial speed [m/s]
th0 = trim_thrust_lo; % trim thrust setting [lb]
de0 = trim_control_lo(1); % trim elevator deflection [rad]

gclear = 40;    % ground clearance [m]

% control saturation limits

ele_lowlim = -25    - de0;
ele_uplim  = 25     - de0;
th_lowlim  = 1000   - th0;
th_uplim   = 19000  - th0;

% -------------------LQR Controller----------------------------------------

% assign weights
w_h  = 100;
w_th = 1;
w_v  = 1;
w_a  = 1;
w_q  = 1;

% assemble Q and R matrix
Q_f = diag([w_h, w_th, w_v, w_a, w_q]);
R_f = diag([0.001 10]);

K_f = lqr(A,B,Q_f,R_f); % full outer matrix

% crop outer LQR matrix
K_o = K_f(:,1);

% crop inner LQR matrix
K_i = K_f(:,2:5);

% -------------------run simulation----------------------------------------

open_system('terrainfollow_sim')
sim('terrainfollow_sim')











