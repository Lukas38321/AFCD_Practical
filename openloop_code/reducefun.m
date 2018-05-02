% load low-fidelity state space matrices created by FindF16Dynamics.m
load sep_lofi_ss

% remove unncessary states
A_lon_OL = A_longitude_lo(2:7,2:7);
B_lon_OL = B_longitude_lo(2:7,:);
