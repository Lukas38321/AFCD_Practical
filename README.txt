Please open the following MATLAB scripts/Simulink models for the corresponding chapters:

Chapter 2: Trim & Linearisation

First the trim and linearisation is done using Matlab script runLINF16sim_ACC.m
Next the accelerometer is added to the model in Simulink model LIN_F16Block_ACC.slx
Analysis of the effect of the accelerometer position is done in \Chapter_5\Chpt5_script.m
	   

Chapter 3: Open Loop Analysis

Chapter 4: Pitch Rate Command System

The 4 state model is reduced to a 2 state model using Matlab script \openloop_code\reducefun_sp.m
The rest of the analysis and design of the pitch-rate control system is done in \openloop_code\openloop_analysis_cpht7.m (note that the first part of this script is simply a copy and paste of the script of Chapter 3 above and the code relevant to this chapter begins on line 227)

Chapter 5: Terrain Following System