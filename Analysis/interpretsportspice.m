% SCRIPT NAME:
%
%   Interpret Tempature Simulation
% 
% DESCRIPTION:
%
%   This script is designed to look at a specific simulation done in
%   LTspice, 'SPORTmodelV5.asc'. The Spice model considers the sport 
%   spacecraft over different spacecraft areas. This data set is used 
%   to verify methods for approximating temperature and density profiles 
%   with function "gettempdensityest.m".
%
% Jason Powell
% Feb 02, 2023

addpath(fullfile('..','utilities'))
boxpath = getboxpath;


%% Load SPORTmodelV5 from .txt file, 
% This only needs to be done once, but is nice to see demonstrate how to get the spice data as a matlab table

sportData = fullfile(boxpath,"LTspice","langmuir","SPORTmodelV5.1.txt");

data = loadspicefile(sportData);
save(fullfile(boxpath,"LTspice","langmuir","SPORTmodelV5.1.mat"),"data")

%% Plot the data

areaarray = unique(data.steps);

for i = 1:height(areaarray)
    d2 = data(data.steps == areaarray(i),:);
    figure
    tiledlayout(2,1)
    % current vs bias voltage
    nexttile
    plot(d2.v1,d2.("Ix(LP:Plasma)"))
    title('Current Voltage Curve')
    xlabel('Bias Voltage [V]')
    ylabel('Current [A]')

    sgtitle(sprintf('SC area: %2.1fm^{2}',d2.steps(i)))

end


%% Attempt to estimate Temperature/ Density

areaarray = unique(data.steps);
density = 1e12;
% for i = 1:height(areaarray)
for i = 3
    d2 = data(data.steps == areaarray(i),:);
    [t,n] = gettempdensityest2(d2.v1,-d2.("Ix(LP:Plasma)"),d2.("V(EFP1,SCground)"),true,1200,density,d2.("V(lp)"));
    
    fprintf('Swenson Temp Guess: %d, Density Guess: %d\n',t(2),n(2))
    fprintf('Swenson Cyl Guess: %d, Density Guess: %d\n',t(3),n(3))
    fprintf('Debchoudhury Temp Guess: %d, Density Guess: %d\n',t(1),n(1))

end