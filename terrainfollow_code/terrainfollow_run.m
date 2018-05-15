%==========================================================================
%
%   Author: Lukas
% 
%==========================================================================
clear;

% load reduced longitudinal state-space matrices and trim state
load redu_ss_terrainfollow % altitude is in m and speed in m/s!
load trim_lo % trim values in ft and ft/s

% simulation parameters

h0  = 1500;  % initial altitude [m]
v0  = 0.3048*trim_state_lo(7); % initial speed [m/s]

th0 = trim_thrust_lo; % trim thrust setting [lb]
de0 = trim_control_lo(1); % trim elevator deflection [deg]

gclear = 40;    % ground clearance [m]

% control saturation limits

ele_lowlim = -25 - de0;
ele_uplim  = 25  - de0;
th_lowlim  = 1000 - th0;
th_uplim  = 19000 - th0;












