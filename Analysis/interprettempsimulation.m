% SCRIPT NAME:
%
%   Interpret Tempature Simulation
% 
% DESCRIPTION:
%
%   This script is designed to look at a specific simulation done in
%   LTspice, 'LangmuirTemp.asc'. The Spice model considers a simple probe
%   model for varying temperature ranges. This data set is used to verify
%   methods for approximating temperature and density profiles with
%   function "gettempdensityest.m".
%
% Jason Powell
% Feb 02, 2023

addpath(fullfile('..','utilities'))
boxpath = getboxpath;


%% Load temperatureData from .txt file, 
% This only needs to be done once, but is nice to see demonstrate how to get the spice data as a matlab table

temperatureData = fullfile(boxpath,"LTspice","langmuir","LangmuirTemp.txt");
data = loadspicefile(temperatureData);
data.appliedpotential = data.v1 + data.("V(scground)");
save(fullfile(boxpath,"LTspice","langmuir","LangmuirTemp.m"),"data")

%% Plot the data

temps = unique(data.steps);

for i = 1:height(temps)
    d2 = data(data.steps == temps(i),:);
    figure
    tiledlayout(2,1)
    % current vs bias voltage
    nexttile
    plot(d2.v1,d2.("Ix(LP:Plasma)"))
    title('Current Voltage Curve')
    xlabel('Bias Voltage [V]')
    ylabel('Current [A]')

    % Applied potential vs bias potential
    nexttile
    plot(d2.v1,d2.v1 + d2.appliedpotential)
    title('Bias Potential vs Applied Potential')
    xlabel('Bias Potential [V]')
    ylabel('Potential Applied to Plasma [V]')

    sgtitle(sprintf('Plasma Temperature: %2.1fK',d2.steps(1)))

end


%% Attempt to estimate Temperature/ Density

temps = unique(data.steps);
density = 1e12;
for i = 1:height(temps)
    d2 = data(data.steps == temps(i),:);
    [t,n] = gettempdensityest(d2.appliedpotential,-d2.("Ix(LP:Plasma)"),true,temps(i),density);
    
    
    fprintf('Swenson Temp Guess: %d, Density Guess: %d\n',t(2),n(2))
    fprintf('Swenson Cyl Guess: %d, Density Guess: %d\n',t(3),n(3))
    fprintf('Debchoudhury Temp Guess: %d, Density Guess: %d\n',t(1),n(1))

end




