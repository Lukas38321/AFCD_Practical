%================================================
%   This program reduces the low-fidelity state-
%   space matrices without actuator dynamics ca-
%   lculated from reducefun.m to the matrix cha-
%   racterising the short period.
%
%   Author: Gervase
% 
%================================================


% load low-fidelity reduced state space matrices created by reducefun.m
load redu_ss

%Velocity assumed to be constant
A_lon_sp = A_lon_ac([3,4],[3,4]);

B_lon_sp = B_lon_ac([3,4],:);

C_lon_sp = C_lon_ac([3,4],[3,4]);

D_lon_sp = D_lon_ac([3,4],:);

save redu_ss_sp A_lon_sp B_lon_sp C_lon_sp D_lon_sp