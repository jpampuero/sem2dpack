function add_subdir_to_path(root_path)
%ADD_SUBDIR_TO_PATH   Add subdirectories to MATLAB path.
%   ADD_SUBDIR_TO_PATH(ROOT_PATH) adds root_path to the MATLAB path
%   variable and then calls itself on all subdirectories of ROOT_PATH.
%
%   Example:
%   add_subdir_to_path(pwd);
%
%   See also SET_PATHS.

%   By Andrew Kositsky
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2009/11/30  $

%%% Add root_path to current MATLAB path
addpath(root_path);

%%% Get list of all items contained in root_path 
files = dir(root_path);

%%% For each item of the list, check if it is a directory and if so call
%%% add_subdir_to_path on it
for i = 1:numel(files)
    if(files(i).isdir == 1) % if the file is a directory...
        if(~ strcmp(files(i).name,'.') && ~ strcmp(files(i).name,'..')) % reject the directories '.' and '..'
            add_subdir_to_path([root_path,'/', files(i).name]); % call this script on the new directory
        end
    end
end