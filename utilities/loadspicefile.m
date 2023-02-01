% Function NAME:
%
%   Load Spice File 
% 
% DESCRIPTION:
%
%   Read in spice data as a .txt file and save the results as a matlab
%   table
%
% INPUTS:
%
%   loc - (string) the absolute path to .txt file of spice data
%
% OUTPUTS:
%
%   data - (table) contains the data read in from text file as matlab table
% 
% ASSUMPTIONS AND LIMITATIONS:
%
%   none
%

%%
% Jason Powell 09 September 2022

function data =  loadspicefile(loc)

    % First read in all the data as strings (helpful for getting step info)
    t = readlines(loc);
    
    % Extract the paramaters from the first line of file
    params = strsplit(t(1));

    % some magicy goo to extract step information
    stepInformation = t(contains(t,"Step Information"));
    
    steps = strings(height(stepInformation),1);
    
    for i = 1:height(stepInformation)
        splits = strsplit(stepInformation(i),{'=',' '});
        steps(i) = splits(4);
    end

    kfactors = contains(steps,"K");
    steps = strrep(steps,'K','');
    steps = str2double(steps);
    steps(kfactors) = steps(kfactors)*1e3;
    
    steps = repelem(steps,501,1);
    
    % Next read in data as a table (gets all the numbers)
    data = readtable(loc);
    % step info cannot be well parsed, so we remove the NaN remnants of them
    data = data(~isnan(data{:,1}),:);

    % rename variables to their parameter names from LTspice
    data.Properties.VariableNames = params;
    % Add the step info
    data = addvars(data,steps,'before',1);
end




















