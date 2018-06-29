%==========================================================================
%   This programm calculates the inherent motion
%   characteristics of the reduced F-16 model
%   - periodic: short period, phugoid, Dutch roll
%   - aperiodic: aperiodic roll, spiral
%  
%   And compares 4-state to 2-state
%
%   Author: Lukas & Gervase
% 
%==========================================================================
clear; 

% import matrices
load redu_ss 
load redu_ss_sp


% define trim point
v0 = 600; %[ft/s]

% set plotting/display parameters
pzmaps          = false;
characteristics = false;
simulations     = false;
comparison      = true;

% ------------------SS matrices--------------------------------------------

% LONGITUDINAL
states_lon    = {'theta' 'v' 'alpha' 'q'};
inputs_lon    = {'elevator'};
outputs_lon   = {'theta' 'v' 'alpha' 'q'};
SS_lon = ss(A_lon, B_lon, C_lon, D_lon,...
            'statename',states_lon,...
            'inputname',inputs_lon,...
            'outputname',outputs_lon);
        
% LATERAL
states_lat    = {'phi' 'beta' 'roll rate' 'yaw rate'};
inputs_lat    = {'aileron' 'rudder'};
outputs_lat   = {'phi' 'beta' 'roll rate' 'yaw rate'};
SS_lat = ss(A_lat, B_lat, C_lat, D_lat,...
            'statename',states_lat,...
            'inputname',inputs_lat,...
            'outputname',outputs_lat);

% -------------------pole-zero maps----------------------------------------
if pzmaps
    
    % LONGITUDINAL
    figure(1); 
    pzmap(SS_lon);
    print -depsc2 -r1200 Figures/pzmap_sym
    
    % LATERAL
    figure(2); 
    pzmap(SS_lat);
    print -depsc2 -r1200 Figures/pzmap_asym
    
end

% ------------------characteristics----------------------------------------
if characteristics
    
    % calculate natural frequencies and damping ratios
    [Wn_lon, zeta_lon, Po_lon] = damp(SS_lon);
    [Wn_lat, zeta_lat, Po_lat] = damp(SS_lat);
    
    % LONGITUDINAL
    
    % calculate period (PDF P.5 check!)
    P_lon = 2*pi()./(Wn_lon.*sqrt(1-zeta_lon.^2));
    % calculate time to half amplitude
    T_halfa_lon = log(2)./(Wn_lon.*zeta_lon);
    % combine characteristics into matrix
    chara_lon = [real(Po_lon),imag(Po_lon),Wn_lon,zeta_lon, P_lon, T_halfa_lon];
    
    % display matrix
    disp('Longitudinal:')
    disp('   Pole (real,imag)       Wn       zeta      P      T_halfa')
    disp('   --------------------------------------------------------')
    disp(chara_lon)
    
    % wrtie matrix to excel
    %xlswrite('Tables/chara_lon',chara_lon);
    
    % LATERAL
    
    % calculate time constant
    tau_lat = 1./(zeta_lat.*Wn_lat);
    % calculate time to half amplitude
    T_halfa_lat = log(2)./(Wn_lat.*zeta_lat);
    % combine characteristics into matrix
    chara_lat = [real(Po_lat), imag(Po_lat), Wn_lat,zeta_lat, tau_lat, T_halfa_lat];
    
    % display matrix
    disp('Lateral:')
    disp('   Pole (real,imag)       Wn       zeta      tau      T_halfa')
    disp('   ----------------------------------------------------------')
    disp(chara_lat)
    
    % wrtie matrix to excel
    %xlswrite('Tables/chara_lat',chara_lat);
    
end

% ---------------------simulation------------------------------------------

if simulations
    
    % figure parameters for printing
    simfig_width    = 1000;
    simfig_height   = 400;
    
    % lsim parameters for pulse-shaped inputs
    T = 10; % observaton period [s]
    dt = 0.01; % sampling time [s]
    
    % pulse parameters
    d_ail = 1; % aileron pulse duration [s]
    m_ail = rad2deg(0.025); % aileron pulse magnitude [deg]
    
    d_rud = 1; % rudder pulse duration [s]
    m_rud = rad2deg(0.025); % rudder pulse magnitude [deg]
    
    % define time axis
    t  = 0:dt:T; N = length(t);
    
    % create signals
    N_ail = d_ail/dt;
    N_rud = d_rud/dt;
    
    sig_ail  = [m_ail*ones(1,N_ail),zeros(1,N-N_ail)];
    sig_rud  = [m_rud*ones(1,N_rud),zeros(1,N-N_rud)];
    sig_null = zeros(1,N);
    
    
    % SHORT PERIOD
    [y_sp, t_sp] = step(-0.05*SS_lon,10); % elevator input: -0.05 [rad]
    theta_sp = y_sp(:,1);
    q_sp = y_sp(:,4);
    
    figure('pos',[100 100 simfig_width simfig_height/2])
    figure(1)
    plot(t_sp,q_sp)
    xlabel('t [sec]')
    ylabel('q [deg/s]')
    grid on
    print -depsc2 -r1200 Figures/short_period
    
    % PHUGOID
    [y_ph, t_ph] = step(-0.05*SS_lon,300); % elevator input: -0.05 [rad]
    v_ph = y_ph(:,2) + v0; % trim velocity added
    alpha_ph = y_ph(:,3);
    
    figure('pos',[100 100 simfig_width simfig_height])
    figure(2)
    subplot(2,1,1)
    plot(t_ph,v_ph)
    xlabel('t [sec]')
    ylabel('V [ft/sec]')
    grid on
    
    subplot(2,1,2)
    plot(t_ph,alpha_ph)
    xlabel('t [sec]')
    ylabel('\alpha [deg]')
    grid on
    print -depsc2 -r1200 Figures/phugoid

    
    % DUTCH ROLL
    % use rudder pulse
    sig_dr = [sig_null',sig_rud'];
    y_dr = lsim(SS_lat,sig_dr,t);
    p_dr = y_dr(:,3);
    r_dr = y_dr(:,4);
    
    figure('pos',[100 100 simfig_width simfig_height])
    figure(3)
    subplot(2,1,1)
    plot(t,p_dr)
    xlabel('t [sec]')
    ylabel('p [deg/sec]')
    grid on
    
    subplot(2,1,2)
    plot(t,r_dr)
    xlabel('t [sec]')
    ylabel('r [deg/s]')
    grid on
    print -depsc2 -r1200 Figures/dutch_roll
    
    % APERIODIC ROLL
    % use aileron pulse
    sig_ap = [sig_ail',sig_null'];
    y_ap = lsim(SS_lat,sig_ap,t);
    p_ap = y_ap(:,3);
    
    % (Plot combined with spiral)
    
    % SPIRAL
    % use aileron pulse
    
    sig_spir = [sig_ail',sig_null'];
    y_spir = lsim(SS_lat,sig_spir,t);
    phi_spir = y_spir(:,1);
    
    figure('pos',[100 100 simfig_width simfig_height])
    figure(4)
    subplot(2,1,1)
    plot(t,p_ap)
    xlabel('t [sec]')
    ylabel('p [deg/sec]')
    grid on
    
    subplot(2,1,2)
    plot(t,phi_spir)
    xlabel('t [sec]')
    ylabel('\phi [deg]')
    grid on
    print -depsc2 -r1200 Figures/aperoll_spiral
    
end

%% Chapter 7
%------------------Short period comparison between 4-state & 2-state-------

if comparison
    
    %2-state SS system
    states_lon_sp    = {'alpha' 'q'};
    inputs_lon_sp    = {'elevator'};
    outputs_lon_sp   = {'alpha' 'q'};
    SS_lon_2s = ss(A_lon_sp, B_lon_sp, C_lon_sp, D_lon_sp,...
                'statename',states_lon_sp,...
                'inputname',inputs_lon_sp,...
                'outputname',outputs_lon_sp);
            
    %4-state
    [y_sp_4s, t_c] = step(SS_lon, 8);
    alpha_sp_4s = y_sp_4s(:,3);
    q_sp_4s = y_sp_4s(:,4);
    
    %2-state
    [y_sp_2s, t_c] = step(SS_lon_2s, 8);
    alpha_sp_2s = y_sp_2s(:,1);
    q_sp_2s = y_sp_2s(:,2);
    
    
  
    %plots
    
    figure(1)
    plot(t_c,q_sp_4s,t_c,q_sp_2s)
    xlabel('t [sec]')
    ylabel('q [deg/s]')
    legend('4-state model', '2-state model', 'Location', 'northeast')
    grid on
    
    
    figure(2)
    subplot(2,1,1)
    plot(t_c,alpha_sp_4s,t_c,alpha_sp_2s)
    legend('4-state model', '2-state model', 'Location', 'northeast')
    xlabel('t [sec]')
    ylabel('\alpha [deg]')
    grid on
    
    subplot(2,1,2)
    plot(t_c,q_sp_4s,t_c,q_sp_2s)
    xlabel('t [sec]')
    ylabel('q [deg/s]')
    legend('4-state model', '2-state model', 'Location', 'northeast')
    grid on
    
end

trans_matrix_redu = tf(SS_lon_2s);
q_to_ele = trans_matrix_redu(2, 1);


[num, den] = tfdata(SS_lon_2s);
K_a = num{1}(3);
K_q = num{2}(3);

T_theta_1 = num{2}(2)/K_q;
w_n_sp_1 = sqrt(den{2}(3));

pol = pole(SS_lon_2s);
[Wn, Zeta] = damp(SS_lon_2s);

W_nsp = 0.03 * v0 * 0.3048; %V in m/s (i.e. v0 * 0.3048)
Z_sp = 0.5;
T_theta_2 = 1/(0.75*W_nsp);

poles = [-Z_sp * W_nsp + W_nsp* sqrt(1-(Z_sp)^2)*i;
        -Z_sp * W_nsp - W_nsp* sqrt(1-(Z_sp)^2)*i];

K = place(A_lon_sp, B_lon_sp, poles);

V_gust = 4.572;          % Design vertical gust, [m/s]
V = v0 * 0.3048;
d_alpha = atan(V_gust/V);

d_del_ele = K(1) * d_alpha;

A_new = A_lon_sp - B_lon_sp * K;

states_lon_sp    = {'alpha' 'q'};
inputs_lon_sp    = {'elevator'};
outputs_lon_sp   = {'alpha' 'q'};
SS_new = ss(A_new, B_lon_sp, C_lon_sp, D_lon_sp,...
            'statename',states_lon_sp,...
            'inputname',inputs_lon_sp,...
            'outputname',outputs_lon_sp);
 
trans_matrix_redu_new = tf(SS_new);
q_to_ele_new = trans_matrix_redu_new(2, 1);

s = tf('s');
F = (T_theta_2*s +1)/(T_theta_1*s + 1); 
K_f = -0.66102;  %obtained from sisotool

final_tf = (7.31874*(s+4.115))/(s^2 + 5.486*s +30.1); %obtained from sisotool

%with prefilter
[y_with_fil, t_fil] = step(final_tf, 2);
%without prefilter
[y_without_fil, t_fil] = step(q_to_ele_new, 2);

figure(3)
plot(t_fil, y_with_fil, t_fil, y_without_fil)
legend('with prefilter', 'without prefilter', 'Location', 'east')
xlabel('t [sec]')
ylabel('q [deg/s]')
grid on


CAP = (W_nsp*W_nsp)/((V/9.80665)*(1/T_theta_2));
DB = T_theta_2 - ((2*Z_sp)/W_nsp);


%ramp response

t_r = 0:0.01:5;
ramp_h = 2;
ramp_in = t_r;
for i = 1:length(t_r)
    if ramp_in(i) >= ramp_h
        ramp_in(i) = ramp_h;
    end
end

model = final_tf;

[y_r,t_r]=lsim(model,ramp_in,t_r);
%[y_ramp,t_r] = step(final_tf/s, 5)

DB_c = max(y_r) - ramp_h;

figure(4)
%plot(t_r, y_ramp)
plot(t_r,y_r, t_r,ramp_in)
hold on
plot(xlim, [1 1]*max(y_r), ':k')
plot(xlim, [1 1]*ramp_h, ':k')
hold off
legend('Pitch angle response', 'ramp input', 'Location', 'southeast')
xlabel('t [sec]')
ylabel('\theta [deg]')
grid on

DB_qs = (max(y_r) - ramp_h)/ramp_h;

%step response
step_h = 1;
step_in = zeros(1, length(t_r));
for i = 1:length(t_r)
    if t_r(i) >= 2
        step_in(i) = step_h;
    end
end

[y_s,t_r]=lsim(model,step_in,t_r);

qm_qss = max(y_s)/step_h;

figure(5)
plot(t_r,y_s, t_r, step_in)
hold on
plot(xlim, [1 1]*max(y_s), ':k')
plot(xlim, [1 1]*step_h, ':k')
hold off
legend('Pitch rate response', 'step input', 'Location', 'southeast')
xlabel('t [sec]')
ylabel('q [deg/s]')
grid on



%% ######### Not sure about this yet ###########
%%

%{
trans_matrix_redu = tf(SS_lon_2s);
q_to_ele = trans_matrix_redu(2, 1);


[num, den] = tfdata(SS_lon_2s);

K_a = num{1}(3);

K_q = num{2}(3);
T_theta2 = num{2}(2)/K_q;
w_n_sp = sqrt(den{2}(3));

pol = pole(SS_lon_2s);
[Wn, Zeta] = damp(SS_lon_2s);

%}

