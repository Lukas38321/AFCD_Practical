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



C_an_1 = [];
for i = 1:length(C_an)
    if C_an(i) ~= 0
        C_an_1 = [C_an_1, C_an(i)];
    end
end
C_an_r = C_an_1(3:end);

%% Question 5: Elevator-to-normal_acceleration transfer function
%%

SS_chpt5 = ss(A_lo, B_lo, C_lo, D_lo);

trans_matrix = tf(SS_chpt5);
ele_to_an_tf = trans_matrix(19,2);


tol = sqrt(eps);
sysr = minreal(SS_chpt5, tol);
trans_matrix = tf(sysr);
ele_to_an_tf_r = trans_matrix(19,2);


%% Question 6: Normal accel. response to step elevator
%%

t1 = 0:0.1:400;
t2 = 0:0.01:5;

y1 = step(-1*ele_to_an_tf, t1);
y2 = step(-1*ele_to_an_tf, t2);


figure(1)
subplot(2,1,1)
plot(t1, y1)
axis([0 400 -0.5 0.65])
xlabel('Time, t [sec]')
ylabel('Normal acceleration, a_{n} [g]')
grid on

subplot(2,1,2)
plot(t2, y2);
axis([0 5 -0.05 0.65])
xlabel('Time, t [sec]')
ylabel('Normal acceleration, a_{n} [g]')
grid on

figure(2)
pzmap(ele_to_an_tf)

%% Question 9: Responses for different accel x-locations
%%

t3 = 0:0.0001:0.2
x_accel = [0, 5, 5.9, 6, 7, 15];
y = [];


for i = 1:length(x_accel)
    x_accel(i)
    legend_info{i} = ['x_{a} = ' num2str(x_accel(i)) ' ft'];
    load(['variables_chpt5_' num2str(x_accel(i)) '.mat'])
    SS_chpt5_total = ss(A_lo, B_lo, C_lo, D_lo);
    trans_matrix_total = tf(SS_chpt5_total);
    ele_to_an_tf_total = trans_matrix_total(19,2);
    [z, gain] =zero(ele_to_an_tf_total)
    pole(ele_to_an_tf_total)
    y = [y, step(-1*ele_to_an_tf_total, t3)];
end



figure(3)
plot(t3, y)
legend(legend_info, 'Location', 'southeast')
xlabel('Time, t [sec]')
ylabel('Normal acceleration, a_{n} [g]')
grid on


figure(4)
hold on;
for i = 1:length(x_accel)
    legend_info_1{i} = ['x_{a} = ' num2str(x_accel(i)) ' ft'];
    load(['variables_chpt5_' num2str(x_accel(i)) '.mat'])
    SS_chpt5_total = ss(A_lo, B_lo, C_lo, D_lo);
    trans_matrix_total = tf(SS_chpt5_total);
    ele_to_an_tf_total = trans_matrix_total(19,2);
    pzmap(ele_to_an_tf_total)
    legend(legend_info_1, 'Location', 'southeast')
    %legend('show')
end
%legend(legend_info_1, 'Location', 'southeast')
hold off;

