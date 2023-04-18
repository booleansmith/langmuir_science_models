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
    
    settings = table;

    % extract setting names
    splits = strsplit(stepInformation(1),{' '});
    settingsList = splits(contains(splits,'='));

    for i = 1:length(settingsList)
        curSetting = strsplit(settingsList(i),'=');
        settings.(curSetting(1)) = strings(height(stepInformation),1);
    end
    
    for i = 1:height(stepInformation)
        splits = strsplit(stepInformation(i),{' '});
        settingsList = splits(contains(splits,'='));
        for j = 1:length(settingsList)
            curSetting = strsplit(settingsList(j),'=');
            settings.(curSetting(1))(i) = curSetting(2);
        end
    end

    % Spice uses SI unit prefixes that matlab cannot interpret, we strip
    % them and add values back later
    for i = 1:width(settings)
        kfactors = contains(settings{:,i},"K");
        Gfactors = contains(settings{:,i},"G");
        mfactors = contains(settings{:,i},"m");
        Tfactors = contains(settings{:,i},"T");
    
        settings{:,i} = strrep(settings{:,i},'K','');
        settings{:,i} = strrep(settings{:,i},'m','');
        settings{:,i} = strrep(settings{:,i},'G','');
        settings{:,i} = strrep(settings{:,i},'T','');
    
        settings.(i) = str2double(settings{:,i});
        settings{kfactors,i} = settings{kfactors,i}*1e3;
        settings{mfactors,i} = settings{mfactors,i}*1e-3;
        settings{Gfactors,i} = settings{Gfactors,i}*1e9;
        settings{Tfactors,i} = settings{Tfactors,i}*1e12;  
        
    end

    settings = repelem(settings,501,1);
    % Next read in data as a table (gets all the numbers)
    data = readtable(loc);
    % step info cannot be well parsed, so we remove the NaN remnants of them
    data = data(~isnan(data{:,1}),:);

    % rename variables to their parameter names from LTspice
    data.Properties.VariableNames = params;
    % Add the step info
    for i = width(settings):-1:1
        data = addvars(data,settings.(i),'before',1,'NewVariableNames',settings.Properties.VariableNames{i});
    end
end




















