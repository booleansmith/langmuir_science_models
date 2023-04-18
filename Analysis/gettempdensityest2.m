function [t,n] = gettempdensityest2(V,I,probe,fpp,dodebug)
% Function NAME:
%
%   Get Temperature and Density Estimation 2
% 
% DESCRIPTION:
%
%   Takes current voltage current measurements and returns temperature and
%   denisty parameters. Version 2 consideres data from fpp's to use in
%   calculating the floating potential of the spacecraft
%
% INPUTS:
%
%   V   - (double array) The bias potential applied to the spacecraft
%   I   - (double array) The measured current for each biased potential
%   fpp - (double array) The floating potential as measured by fpp
%
% OPTIONAL INPUTS
%
%   dodebug - (boolean) optionally enable plots for debugging
%   temp    - (double) the actual temperature
%   density - (double) the actual density
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
% Feb 02, 2023
%%

arguments
    V (1,:) double
    I (1,:) double
    probe (1,1) struct
    fpp (1,:) double = zeros(size(V))
    dodebug logical = false;
end



% KNOWN CONSTANTS
k_b = 1.38e-23;  % Boltsmann constant m^2*kg*s^-2*K^-1
e   = 1.602e-19; % elementary charge
m_e = 9.109e-31; % electron mass, kg
m_i = 2.65e-26; % O+ ion mass, kg


% Probe/ SC parameters
probeRadius = 0.5e-3;    % meters
probeHeight = 85.2e-3; % meters
A_proj = probeRadius*2*probeHeight; % projected area, m^2
% A      = 2*pi*probeRadius*probeHeight; % OML collection area, m^2
A = 297.287e-6; % TODO: The actual area of SPORT LP needs to be decided
beta = 0.5;
V_sc  = 7.4e3; % SC velocity, m/s



%% Pre-treat the data a little bit
% This is necessary to interpolate values out of IV curve, hopefully no
% meaningful data is lost

[I,ia,~] = unique(I,'stable');
V = V(ia);

%% Solve for stage 1 parameters

% floating potential
vf = interp1(I,V,0); % The voltage at the zero current threshold

% plasma potential
Idot = gradient(I,V); % di/dv
Idot_max = max(Idot(V > vf & V < vf+1.5));
phi_p = V(Idot == Idot_max);

% ion density
index = V > -2 & V < -1.5;
p = polyfit(V(index),I(index),1);
I_ram_est = -polyval(p,phi_p);
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
    plot(phi_p,Idot_max,"*",'MarkerSize',5)
    xlabel('Bias Potential (V)')
    ylabel('dI/dV')

    title('First Stage parameters from IV curve')
    legend('I','\Phi_f','','','I_{Ram_{i}}','dI/dV','\Phi_p')
end

%% Estimate the sheath potential

% The floating potential meters monitor the change in the SC potential,
% this change needs to be included in sheath potential considerations
minFpp = min(fpp);
deltaVf = fpp - minFpp;

% Calculate the potential across the sheath
phi_s = V - phi_p - fpp;
% phi_s = V - fpp;


% if dodebug
%     figure
%     plot(V,appliedPotential)
%     hold on
%     plot(V,phi_s)
%     title('Verification of Sheath Potential Estimation')
%     legend('Actual','Estimation')
%     xlabel('Probe Potential, \phi_0')
%     ylabel('Sheath Potential, \phi_s')
% end


%% Do non-linear curve fitting

% x(1): n - number density
% x(2): T - temperature
% x(3): A - Area of the probe
% x(4): Rp - Radius of the probe
% x(5): w     - speed of the probe relative to the plasma (m/s)
% x(6): theta - the angle between the probe axis and the velocity
%               vector (radians)
% x(7): beta 

fun = @(x,phi) -x(1)*A_proj*e*V_sc ... % ion ram current
               -x(1)*e*A*sqrt((k_b*x(2))/(2*pi*m_i))*(1-((e*phi)/(k_b*x(2)))).^x(7) ... % ion current
               -x(1)*-e*A*sqrt((k_b*x(2))/(2*pi*m_e))*exp((e*phi)/(k_b*x(2))); % electron current
% Ion saturation region/ Electron retradation region mask
m = phi_s < 0.01; 
m2 = appliedPotential < 0.01;


T_guess = 1250; % kelvin
T_max = 5000;
T_min = 300;

x0 = [ni_est,T_guess,A,probeRadius, V_sc, pi/2,beta]; % initial guesses for lsqcurvefit
bounds(1,:) = [ni_est - 0.5*ni_est, T_min, A, probeRadius, V_sc, pi/2,0];
bounds(2,:) = [ni_est + 0.5*ni_est, T_max, A, probeRadius, V_sc, pi/2,1];
tolerance = 1e-25;
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','StepTolerance',tolerance,'FunctionTolerance',tolerance);
% options = optimoptions('lsqcurvefit','Display','none');
x1 = lsqcurvefit(fun                 ,x0,phi_s(m),I(m),bounds(1,:),bounds(2,:),options);
x2 = lsqcurvefit(@OMLCurrentApprox   ,x0,phi_s,I,bounds(1,:),bounds(2,:),options);
x3 = lsqcurvefit(@OMLCurrentCyl      ,x0,phi_s,I,bounds(1,:),bounds(2,:),options);
x4 = lsqcurvefit(@OMLCurrentCylMoving,x0,phi_s,I,bounds(1,:),bounds(2,:),options);

xtrue = [density,temp,A,probeRadius, V_sc, pi/2, beta];


if dodebug
    
    figure
    tiledlayout('flow')

    nexttile % Debchoudhury
    semilogy(phi_s,abs(I))
    hold on
    semilogy(phi_s(m),abs(fun(x1,phi_s(m))))
    semilogy(phi_s(m),abs(fun(xtrue,phi_s(m))))
    legend('LTSpice','Line Fit','Model','Location','southwest')
    xlabel('Potential Across Sheath, \phi_s [V]')
    ylabel('Current Collected by Probe [A]')
    title('Debchoudhury')
    str = sprintf('N_{est} = %0.2e, T_{est} = %0.2e', x1(1),x1(2));
    text(min(phi_s)*9/10,0.5*mean(I),str);

    nexttile % Swenson OML Approx
    semilogy(phi_s,abs(I))
    hold on
    semilogy(phi_s,abs(OMLCurrentApprox(x2,phi_s)))
    semilogy(phi_s,abs(OMLCurrentApprox(xtrue,phi_s)))
    legend('LTSpice','Line Fit','Model','Location','southwest')
    xlabel('Potential Across Sheath, \phi_s [V]')
    ylabel('Current Collected by Probe [A]')
    title('Swenson OML Approx')
    str = sprintf('N_{est} = %0.2e, T_{est} = %0.2e', x2(1),x2(2));
    text(min(phi_s)*9/10,0.5*mean(I),str);
    
    nexttile % Swenson OML
    semilogy(phi_s,abs(I))
    hold on
    semilogy(phi_s,abs(OMLCurrentCyl(x3,phi_s)))
    semilogy(phi_s,abs(OMLCurrentCyl(xtrue,phi_s)))
    legend('LTSpice','Line Fit','Model','Location','southwest')
    xlabel('Potential Across Sheath, \phi_s [V]')
    ylabel('Current Collected by Probe [A]')
    title('Swenson OML Cyl')
    str = sprintf('N_{est} = %0.2e, T_{est} = %0.2e', x3(1),x3(2));
    text(min(phi_s)*9/10,0.5*mean(I),str);

    nexttile % Swenson OML Moving
    semilogy(phi_s,abs(I))
    hold on
    semilogy(phi_s,abs(OMLCurrentCylMoving(x4,phi_s)))
    semilogy(phi_s,abs(OMLCurrentCylMoving(xtrue,phi_s)))
    legend('LTSpice','Line Fit','Model','Location','southwest')
    xlabel('Potential Across Sheath, \phi_s [V]')
    ylabel('Current Collected by Probe [A]')
    title('Swenson OML Cyl Moving')
    str = sprintf('N_{est} = %0.2e, T_{est} = %0.2e', x3(1),x3(2));
    text(min(phi_s)*9/10,0.5*mean(I),str);

    titlestring = sprintf('Actual Paramters: n = %0.2e, T = %0.2e',density,temp);
    sgtitle(titlestring)
    

end
t = [x1(2),x2(2),x3(2)];
n = [x1(1),x2(1),x3(1)];

end

