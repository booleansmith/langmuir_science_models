% PROCESS NAME:
%
%   Langmuir Probe Data
% 
% DESCRIPTION:
%
%   Produce syntetic Langmuir probe data to analyse
%
% INPUTS:
%
%   none
%
% OUTPUTS:
%
%   data - (table) contains entries for LP current measurements and 
%          respective voltage potential levels
% 
% ASSUMPTIONS AND LIMITATIONS:
%
%   none
%


% Jason Powell 09 September 2022

% The first model we consider is based on Debchoudhury's method for looking
% at FPMU data on the ISS: "Observations and Validation of Plasma Density,
% Temperature, and O+ Abundance From a Langmuir Probe Onboard the ISS"

% I = -I_ram_i - I_OML_i + I_e
%   I_ram_i = n_j*A_proj*e*V_sc ; The ions that run into collector body
%   I_OML_i = I_th_j*({1-e(phi-phi_p)}/k_b*T_i)^beta ; the OML ion current
%   I_e     = I_th_e*exp({e(phi-phi_p)}/k_b*T_e) ; the electron current 
%   I_th_j  = n_j*q_j*A*sqrt(k_b*T_j/(2*pi*m_j)) ; the thermal current for ions/electrons

% KNOWN CONSTANTS
k_b = 1.38e-23;  % Boltsmann constant m^2*kg*s^-2*K^-1
e   = 1.602e-19; % elementary charge
q_e   = e; % electron charge coulombs 
q_i   = e;  % ion charge coulombs
m_e = 9.109e-31; % electron mass, kg

% Environment Parameters
n = 1e12; % density m^-3
T = 1160; % Temperature, K
m_i = 2.65e-26; % O+ ion mass, kg

phi_p = 1; % plasma potential, Volts
phi = -2:0.01:phi_p; % Bias potential, Volts
V_sc  = 7.4e3; % SC velocity, m/s


% Probe parameters
% probeRadius = 0.5e-3;    % meters
% probeHeight = 85.2e-3; % meters
probeRadius = 0.5e-3;    % meters
probeHeight = 85.2e-3; % meters
A_proj = probeRadius*2*probeHeight; % projected area, m^2
A      = pi*probeRadius*probeHeight; % OML collection area, m^2
beta = 1;

I_ram  = n*A_proj*e*V_sc;
I_th_e = n*e*A*sqrt((k_b*T)/(2*pi*m_e));
I_th_i = n*e*A*sqrt((k_b*T)/(2*pi*m_i));
I_OML_i = I_th_i*(1-((e*(phi-phi_p))/(k_b*T))).^beta;
I_e     = I_th_e*exp((e*(phi-phi_p))/(k_b*T));

I = -I_ram -I_OML_i + I_e;

figure
plot(phi,I)
xlabel('Probe Potential [V]')
ylabel('Current out of probe [A]')
title('Langmuir Probe IV curve')

%% Solve for stage 1 parameters

% floating potential
vf = interp1(I,phi,0); % The voltage at the zero current threshold

% fit a line to the ion current collection
index = phi > -2 & phi < -1.5;
p = polyfit(phi(index),I(index),1);
I_ram_est = -polyval(p,phi_p);
n_est = I_ram_est/(A_proj*e*V_sc);

figure
plot(phi,I)
hold on
plot(phi,polyval(p,phi))
plot(phi_p,-I_ram_est,'r*')
yline(-I_ram_est)
xline(phi_p)
text(phi_p+0.05,-I_ram_est+0.05,['I_{ram} = ', sprintf('%0.2f',I_ram_est)])
xlim([-2,2])
title('Ram Current Estimation')

%% Do non-linear curve fitting

% x(1): n - number density
% x(2): T - temperature
% x(3): phi_p - plasma potential
% x(4): beta - OML parameter
fun = @(x,phi) -x(1)*A_proj*e*V_sc ... % ion ram current
                        -x(1)*e*A*sqrt((k_b*x(2))/(2*pi*m_i))*(1-((e*(phi-x(3)))/(k_b*x(2)))).^x(4) ... % ion current
                        +x(1)*e*A*sqrt((k_b*x(2))/(2*pi*m_e))*exp((e*(phi-x(3)))/(k_b*x(2))); % electron current


% for now we use known values for phi_p and beta
T_guess = 1200; % kelvin
phi_guess = phi_p;
beta_guess = beta;

x0 = [n_est,T_guess,phi_guess,beta_guess]; % initial guesses for lsqcurvefit
bounds(1,:) = x0 - 0.05*x0;
bounds(2,:) = x0 + 0.05*x0;
tolerance = 1e-25;
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','StepTolerance',tolerance,'FunctionTolerance',tolerance);
x = lsqcurvefit(fun,x0,phi,I,bounds(1,:),bounds(2,:),options);


figure
tiledlayout(4,1)
nexttile([3,1])
plot(phi,I)
hold on
plot(phi,fun(x,phi))
nexttile
plot(phi,abs(I-fun(x,phi)))























