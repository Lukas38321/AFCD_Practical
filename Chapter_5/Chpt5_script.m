%================================================
% Chapter 5 script 

%
% Author: Gervase
% 
%================================================

load variables_chpt5_0

%% Question 3: Output equation for a_n
%%

C_an = C_lo(19,:);
D_an = D_lo(19,:);

%% Question 5: Elevator-to-normal_acceleration transfer function
%%

SS_chpt5 = ss(A_lo, B_lo, C_lo, D_lo);

trans_matrix = tf(SS_chpt5);
ele_to_an_tf = trans_matrix(19,2);

%% Question 6: Normal accel. response to step elevator
%%

t = 0:0.01:5;
 
y = step(-1*ele_to_an_tf, t);

%{
figure(1)
subplot(2,1,1)
step(-1*ele_to_an_tf)

subplot(2,1,2)
plot(t, y);
grid on
%}
%% Question 9: Responses for different accel x-locations
%%

x_accel = [0, 5, 5.9, 6, 7, 15];
y = [];
legendInfo = [];

for ii = 1:length(x_accel)
    load(['variables_chpt5_' num2str(x_accel(ii)) '.mat'])
    SS_chpt5_total = ss(A_lo, B_lo, C_lo, D_lo);
    trans_matrix_total = tf(SS_chpt5_total);
    ele_to_an_tf_total = trans_matrix_total(19,2);
    
    y = [y, step(-1*ele_to_an_tf_total, t)];
    
    legendInfo = [legendInfo, ['x_{a} = ' num2str(x_accel(ii)) ' ft,']];
end
legendinfo = strsplit(legendInfo,',');
legendinfo(end) =[];

figure(3)
plot(t, y)
legend(legendinfo, 'Location', 'southeast')
xlabel('Time, t [sec]')
ylabel('Normal acceleration, a_{n} [g]')
grid on

