function [t,n] = gettempdensityest(V,I,dodebug,temp,density)
% Function NAME:
%
%   Get Temperature and Density Estimation
% 
% DESCRIPTION:
%
%   Takes current voltage current measurements and returns temperature and
%   denisty parameters
%
% INPUTS:
%
%   V  - (double array) a vector of applied potentials
%   I  - (double array) a vector of measured currents
%        corresponding to applied potentials
%
% OPTIONAL INPUTS
%
%   dodebug - (boolean) optionally enable plots for debugging
%   temp    - (double) the actual temperature measurement
%   density - (double) the actual density measurement
%
% OUTPUTS:
%
%   t - (double) the temperature estimate produced by function
%   n - (double) the density estimate produced by the function
% 
% ASSUMPTIONS AND LIMITATIONS:
%
%   An oxygen ion plasma is assumed
%
% Jason Powell
% January 31, 2023
%%

if ~exist('dodebug','var')
    dodebug = false;
end

if ~exist('temp','var')
    temp = 1200;
    warning('Temperature data assumed')
end

if ~exist('density','var')
    density = 1e12;
    warning('Density data assumed')
end


% KNOWN CONSTANTS
k_b = 1.38e-23;  % Boltsmann constant m^2*kg*s^-2*K^-1
e   = 1.602e-19; % elementary charge
q_e   = e; % electron charge coulombs 
q_i   = e;  % ion charge coulombs
m_e = 9.109e-31; % electron mass, kg

% Probe parameters
probeRadius = 0.5e-3;    % meters
probeHeight = 85.2e-3; % meters
A_proj = probeRadius*2*probeHeight; % projected area, m^2
% A      = 2*pi*probeRadius*probeHeight; % OML collection area, m^2
A = 297.287e-6; % TODO: The actual area of SPORT LP needs to be decided
Rp = 0.5e-3; % radius of the probe, m
beta = 0.5;
V_sc  = 7.4e3; % SC velocity, m/s
m_i = 2.65e-26; % O+ ion mass, kg


%% Pre-treat the data a little bit

[I,ia,~] = unique(I,'stable');
V = V(ia);

%% Solve for stage 1 parameters

% floating potential
vf = interp1(I,V,0); % The voltage at the zero current threshold

% plasma potential
Idot = gradient(I,V); % di/dv
phi_p = max(Idot(V > vf & V < vf+1.5));

% ion density
index = V > -2 & V < -1.5;
p = polyfit(V(index),I(index),1);
I_ram_est = -polyval(p,phi_p);
% Note, when the sc does not have a velocity, this method will not work.
V_sc = 3.0258e+03; 
ni_est = I_ram_est/(A_proj*e*V_sc);


if dodebug 
    figure
    plot(V,I)
    hold on
    plot(vf,0,"*",'MarkerSize',5)
    plot(V,polyval(p,V),'--k','LineWidth',1)
    xline(phi_p,'--k')
    plot(phi_p,-I_ram_est,"*",'MarkerSize',5)
    text(phi_p+.05,-I_ram_est-2.5e-7,sprintf('n_i = %0.2e',ni_est))
    ylabel('Return Current (A)')
    yyaxis right
    plot(V,Idot)
    hold on
    plot(phi_p,Idot(Idot == phi_p),"*",'MarkerSize',5)
    xlabel('Bias Potential (V)')
    ylabel('dI/dV')

    title('First Stage parameters from IV curve')
    legend('I','Vf','','','I_{Ram}','dI/dV','\Phi_p')
end


%% Do non-linear curve fitting

% x(1): n - number density
% x(2): T - temperature
% x(3): phi_p - plasma potential
% x(4): beta - OML parameter
% x(5): A - Area of the probe
% x(6): Rp - Radius of the probe
fun = @(x,phi) -x(1)*A_proj*e*V_sc ... % ion ram current
               -x(1)*e*A*sqrt((k_b*x(2))/(2*pi*m_i))*(1-((e*phi)/(k_b*x(2)))).^x(3) ... % ion current
               -x(1)*-e*A*sqrt((k_b*x(2))/(2*pi*m_e))*exp((e*phi)/(k_b*x(2))); % electron current
 



% for now we use known values for phi_p and beta
T_guess = 1250; % kelvin
phi_guess = phi_p;
% beta_guess = beta;
beta = 0.5;

% Ion saturation region/ Electron retradation region mask
m = V < 0.01;
x0 = [ni_est,T_guess,beta,A,Rp]; % initial guesses for lsqcurvefit
bounds(1,:) = [ni_est - 0.5*ni_est, 300,phi_guess,beta, A, Rp];
bounds(2,:) = [ni_est + 0.5*ni_est, 5000,phi_guess,beta, A, Rp];
tolerance = 1e-25;
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','StepTolerance',tolerance,'FunctionTolerance',tolerance);
x1 = lsqcurvefit(fun,x0,V(m),I(m),bounds(1,:),bounds(2,:),options);
x2 = lsqcurvefit(@OMLCurrentApprox,x0,V,I,bounds(1,:),bounds(2,:),options);
x3 = lsqcurvefit(@OMLCurrentCyl,x0,V,I,bounds(1,:),bounds(2,:),options);

xtrue = [density,temp,phi_guess,beta,A,Rp];


if dodebug
    
    

    figure
    tiledlayout('flow')

    nexttile % Debchoudhury
    plot(V,I)
    hold on
    plot(V(m),fun(x1,V(m)))
    plot(V(m),fun(xtrue,V(m)))
    legend('LTSpice','Line Fit','Matlab','Location','northwest')
    xlabel('Potential Across Sheath, \phi_s [V]')
    ylabel('Current Collected by Probe [A]')
    title('Debchoudhury')
    str = sprintf('N_{est} = %0.2e, T_{est} = %0.2e', x1(1),x1(2));
    text(0.75,0.5*mean(I),str);

    nexttile % Swenson OML Approx
    plot(V,I)
    hold on
    plot(V,OMLCurrentApprox(x2,V))
    plot(V,OMLCurrentApprox(xtrue,V))
    legend('LTSpice','Line Fit','Matlab','Location','northwest')
    xlabel('Potential Across Sheath, \phi_s [V]')
    ylabel('Current Collected by Probe [A]')
    title('Swenson OML Approx')
    str = sprintf('N_{est} = %0.2e, T_{est} = %0.2e', x2(1),x2(2));
    text(0.75,0.5*mean(I),str);
    
    nexttile % Swenson OML
    plot(V,I)
    hold on
    plot(V,OMLCurrentCyl(x3,V))
    plot(V,OMLCurrentCyl(xtrue,V))
    legend('LTSpice','Line Fit','Matlab','Location','northwest')
    xlabel('Potential Across Sheath, \phi_s [V]')
    ylabel('Current Collected by Probe [A]')
    title('Swenson OML Cyl')
    str = sprintf('N_{est} = %0.2e, T_{est} = %0.2e', x3(1),x3(2));
    text(0.75,0.5*mean(I),str);

    titlestring = sprintf('Actual Paramters: n = %0.2e, T = %0.2e',density,temp);
    sgtitle(titlestring)
    

end
t = [x1(2),x2(2),x3(2)];
n = [x1(1),x2(1),x3(1)];

end

