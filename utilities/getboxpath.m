function [boxPath] = getboxpath()
% FUNCTION NAME:
%
%   getboxpaths
% 
% DESCRIPTION:
%
%   Allows user to specify a path to shared box data. Scripts beyond this 
%   should use this path in to reference data within the box data
%   directory, instead of local paths from their computer. This will allow
%   other users to use the same code without changing the references.
%   Function has an option to remember the box path, so a user does not
%   have to specify this every time. The path is save locally as a .mat
%   file and is not tracked with git.
%
% INPUTS:
%
%   None - 
%
% OUTPUTS:
%
%   boxPath - (char array) Containing the path to box file system
% 
% ASSUMPTIONS AND LIMITATIONS:
%
% None

% Jason Powell July 21st 2022

    % Detect if user has saved box path (in function library)
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    if isfile(fullfile(filepath,'boxpath.mat'))
        load('boxpath.mat','boxPath')
    else
        boxPath = setboxdir;
    end

end

% Function to prompt user to select box directory, and optionally save as a
% variable
function boxPath = setboxdir

    if ispc
        boxPath = uigetdir('C:\','Select the path the the box directory');
    else
        boxPath = uigetdir('~','Select the path the the box directory');
    end

    answer = questdlg('Would you like Matlab to remember your preference?', ...
             'Remember Box Path', ...
             'Yes','No','Yes');

    switch answer
        case 'Yes'
            % Find the function library path (where the variable is stored)
            [filepath,~,~] = fileparts(mfilename('fullpath'));
            save(fullfile(filepath,'boxpath.mat'),'boxPath');
        case 'No'
            return
        otherwise
            return
    end

end