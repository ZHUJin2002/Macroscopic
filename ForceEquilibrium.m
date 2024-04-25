clc
clear
close all

% Parameters
eta = 0.001; % Dynamic viscosity of the fluid (Pa.s)
a = 2.5e-6/2; % Semi-major axis of the ellipsoid (m)
b = 1e-6/2; % Semi-minor axis of the ellipsoid (m)
V_Magnitude = 40e-6; % Pre-set velocity (not the final)

F_preset_parallel = 1e-12;

Proportion = [1/35, 1/5, 1, 4, 68]; %F2/F1
Angle1 = [16.1, 15.0, 16.6, 27.9, 33.3];
Angle2 = [-33.1, -27.2, -16.2, -16.2, -13.6];

F1_parallel = F_preset_parallel./ (Proportion+1);
F1_perp = F1_parallel.*tand(Angle1);

F2_parallel = F_preset_parallel - F1_parallel;
F2_perp = F2_parallel.*tand(Angle2);

F_parallel_total = F1_parallel + F2_parallel;
F_perp_total = F2_perp + F1_perp;

% some basic 
tau0 = 1/sqrt(1 - (b/a)^2);
lambda0 = 1/sqrt((a/b)^2 - 1);

K_parallel = 4/(3*sqrt(tau0^2-1)*((tau0^2+1)*atanh(1/tau0)-tau0));
K_perp = 4/(3*sqrt(lambda0^2+1)*(lambda0-(lambda0^2-1)*acot(lambda0)));


cosofValueAngle = F_parallel_total./(-6*pi*eta*b*K_parallel*V_Magnitude); sinofValueAngle = F_perp_total./(-6*pi*eta*b*K_perp*V_Magnitude);

WobbleAngle = 2*atand(sinofValueAngle./cosofValueAngle)

%% normalization & compute V
actualProportion = sqrt(cosofValueAngle.^2 + sinofValueAngle.^2);
V_norm = V_Magnitude/actualProportion(3); %take 0 phase difference case as norm
F_norm_parallel = -6*pi*eta*b*K_parallel*V_norm*cosd(WobbleAngle(3)/2)

F1_parallel_update = F_norm_parallel./ (Proportion+1);F1_perp_update = F1_parallel.*tand(Angle1);

F2_parallel_update = F_norm_parallel - F1_parallel; F2_perp_update = F2_parallel.*tand(Angle2);

F_parallel_total_update = F1_parallel_update + F2_parallel_update;
F_perp_total_update = F2_perp_update + F1_perp_update;
cosofValueAngle_update = F_parallel_total_update./(-6*pi*eta*b*K_parallel*V_norm); sinofValueAngle_update = F_perp_total_update./(-6*pi*eta*b*K_perp*V_norm);
normLength = sqrt(cosofValueAngle.^2+sinofValueAngle_update.^2);
V_update = V_norm*normLength




