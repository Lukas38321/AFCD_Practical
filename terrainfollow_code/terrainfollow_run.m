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
R_f = diag([0.001 15]);

K_f = lqr(A,B,Q_f,R_f); % full outer matrix

% crop outer LQR matrix
K_o = K_f(:,1);

% crop inner LQR matrix
K_i = K_f(:,2:5);

% -------------------run simulation----------------------------------------

open_system('terrainfollow_sim')
sim('terrainfollow_sim')

% -------------------plot figures------------------------------------------

% altitude figure
figure('pos',[100 100 1200 800]) 

% flightpath(distance)
subplot(2,1,1)
hold on
plot(distance.data(:,1),flightpath.data(:,3),distance.data(:,1),flightpath.data(:,2),':k')
canyon = area(distance.data(:,1),flightpath.data(:,1));
hold off
canyon.FaceColor = [0.5 0.5 0.5];
axis([0,5000,1300,1550])
title('Flight Path')
xlabel('Position [m]')
ylabel('Altitude [m]')
legend('Flight Path','Reference Altitude','Canyon Profile','Location','southeast')
grid on

% altitude error(distance)
subplot(2,1,2)
hold on
alterr = plot(distance.data(:,1),alt_error.data(:,3));
uplim  = plot(distance.data(:,1),alt_error.data(:,2),':k');
lowlim = plot(distance.data(:,1),alt_error.data(:,1),':k');
hold off
axis([0,5000,-20,20])
title('Tracking Performance')
xlabel('Position [m]')
ylabel('Error [m]')
legend([alterr,uplim],{'Altitude Error','Overshoot Limits'},'Location','southeast')
grid on

print -depsc2 -r1200 figures/flightpath_err


% Actuator Figure
figure('pos',[100 100 1200 800])

% thrust setting(time)
subplot(2,1,1)
hold on
thrustset   = plot(thrust.time,thrust.data(:,1));
uplim_th    = plot(thrust.time,thrust.data(:,2),':k');
lowlim_th   = plot(thrust.time,thrust.data(:,3),':k');
hold off
axis([0,60,-2500,17000])
title('Thrust Input')
xlabel('Time [s]')
ylabel('Thrust Setting [lb]')
legend([thrustset,uplim_th],{'Thrust Setting','Saturation Limits'},'Location','northeast')
grid on

% elevator deflection(time)
subplot(2,1,2)
hold on
eleset      = plot(elevator.time,elevator.data(:,1));
uplim_ele   = plot(elevator.time,elevator.data(:,2),':k');
lowlim_ele  = plot(elevator.time,elevator.data(:,3),':k');
hold off
axis([0,60,-30,35])
title('Elevator Input')
xlabel('Time [s]')
ylabel('Deflection [deg]')
legend([eleset,uplim_ele],{'Commanded Deflection','Saturation Limits'},'Location','northeast')
grid on

print -depsc2 -r1200 figures/control_inputs












