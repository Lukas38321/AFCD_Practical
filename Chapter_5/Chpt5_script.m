%================================================
% Chapter 5 script 

%
% Author: Gervase
% 
%================================================

load variables_chpt5

%% Question 3: Output equation for a_n
%%

C_an = C_lo(19,:);
D_an = D_lo(19,:);

%% Question 5: elevator-to-normal_acceleration transfer function
%%

SS_chpt5 = ss(A_lo, B_lo, C_lo, D_lo);


