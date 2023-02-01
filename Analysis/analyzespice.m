% SCRIPT NAME:
%
%   Analyze Spice
% 
% DESCRIPTION:
%
%   Look at LP data computed in SPICE and try to come up with temperature
%   and density parmaters
% 
% ASSUMPTIONS AND LIMITATIONS:
%

% Jason Powell 
% October 2022

%% Set up paths
addpath(fullfile('..','utilities/'))
boxPath = getboxpath();
twobodydata = fullfile(boxPath,'LTspice','Matlab','Test_TwoFloatingBodies_Sweep.mat');

% Load Data
load(twobodydata);

%% Known Parameters

% Probe parameters
% probeRadius = 0.5e-3;    % meters
% probeHeight = 85.2e-3; % meters
probeRadius = 0.5e-3;    % meters
probeHeight = 85.2e-3; % meters
A_proj = probeRadius*2*probeHeight; % projected area, m^2
A      = pi*probeRadius*probeHeight; % OML collection area, m^2
beta = 1;

%% For now, just conisder the largest SC area data

SCarea = 0.062;

S.Table = S.Table(S.Table.scarea == SCarea,:);

plotspicedata(S)

%% Fit the SPICE data to the model

% Function to fit LP data too
fun = @(x,phi) -x(1)*A_proj*e*V_sc ... % ion ram current
                        -x(1)*e*A*sqrt((k_b*x(2))/(2*pi*m_i))*(1-((e*(phi-x(3)))/(k_b*x(2)))).^x(4) ... % ion current
                        +x(1)*e*A*sqrt((k_b*x(2))/(2*pi*m_e))*exp((e*(phi-x(3)))/(k_b*x(2))); % electron current



% Initial Guess
T_guess = 1200; % kelvin
phi_guess = phi_p;
beta_guess = beta;

x0 = [n_est,T_guess,phi_guess,beta_guess]; % initial guesses for lsqcurvefit
bounds(1,:) = x0 - 0.05*x0;
bounds(2,:) = x0 + 0.05*x0;
tolerance = 1e-25;
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','StepTolerance',tolerance,'FunctionTolerance',tolerance);
x = lsqcurvefit(fun,x0,phi,I,bounds(1,:),bounds(2,:),options);