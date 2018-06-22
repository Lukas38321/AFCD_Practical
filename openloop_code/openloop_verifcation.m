%==========================================================================
%   This programm verifies the matrices created by reducefun_openloop by
%   comparing the input/output relation to the original
%
%   Author: Lukas
% 
%==========================================================================
clear;

% import state-space matrices
reducefun_openloop

% run simulation
sim('reducefun_verification')

% create figures

% Longitudinal: Elevator
theta_r   = elevator_out.data(:,1); 
v_r       = elevator_out.data(:,2); 
alpha_r   = elevator_out.data(:,3);  
q_r       = elevator_out.data(:,4);
theta_o   = elevator_out.data(:,5);
v_o       = elevator_out.data(:,6);
alpha_o   = elevator_out.data(:,7); 
q_o       = elevator_out.data(:,8); 

figure('pos',[100 100 1200 400])
plot(tout,theta_r, tout,theta_o,'x',...
    tout,v_r,   tout,v_o,'x',...
    tout,alpha_r, tout,alpha_o,'x',...
    tout,q_r,tout,q_o,'x');
legend('theta r [deg]','theta o [deg]',...
    'vt r [ft/s]','vt o [ft/s]',...
    'alpha r [deg]','alpha o [deg]',...
    'q r [deg/s]','q o [deg/s]',...
    'Location','southeast')
xlabel('time [s]')
print -depsc2 -r1200 figures/veri_elevator


% Lateral: Aileron
phi_r   = aileron_out.data(:,1); 
beta_r  = aileron_out.data(:,2); 
p_r     = aileron_out.data(:,3); 
r_r     = aileron_out.data(:,4);
phi_o   = aileron_out.data(:,5); 
beta_o  = aileron_out.data(:,6); 
p_o     = aileron_out.data(:,7); 
r_o     = aileron_out.data(:,8);

figure('pos',[100 100 1200 400])
plot(tout,phi_r,    tout,phi_o,'x',...
    tout,beta_r,    tout,beta_o,'x',...
    tout,p_r,       tout,p_o,'x',...
    tout,r_r,       tout,r_o,'x');
legend('phi r [deg]', 'phi o [deg]',...
       'beta r [ft/s]','beta o [ft/s]',...
       'p r [deg]','p o [deg]',...
       'r r [deg/s]','r o [deg/s]',...
       'Location','southeast')
xlabel('time [s]')
print -depsc2 -r1200 figures/veri_aileron


% Lateral: Rudder
phi_r   = rudder_out.data(:,1); 
beta_r  = rudder_out.data(:,2); 
p_r     = rudder_out.data(:,3); 
r_r     = rudder_out.data(:,4);
phi_o   = rudder_out.data(:,5); 
beta_o  = rudder_out.data(:,6); 
p_o     = rudder_out.data(:,7); 
r_o     = rudder_out.data(:,8);

figure('pos',[100 100 1200 400])
plot(tout,phi_r,    tout,phi_o,'x',...
    tout,beta_r,    tout,beta_o,'x',...
    tout,p_r,       tout,p_o,'x',...
    tout,r_r,       tout,r_o,'x');
legend('phi r [deg]', 'phi o [deg]',...
       'beta r [ft/s]','beta o [ft/s]',...
       'p r [deg]','p o [deg]',...
       'r r [deg/s]','r o [deg/s]',...
       'Location','southeast')
xlabel('time [s]')
print -depsc2 -r1200 figures/veri_rudder

