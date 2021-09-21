function set_paths(code_dir)
%SET_PATHS   Set MATLAB path for running PCAIM_driver.
%   SET_PATHS(CODE_DIR) adds the Code directory of the PCAIM software
%   package to the current MATLAB path. Note that this does NOT change your
%   default path settings and the command must be rerun each time MATLAB is
%   opened.
%   
%   Example:
%   code_dir = pwd;
%   set_paths(code_dir);
%
%   See also ADD_SUBDIR_TO_PATH, PCAIM_DRIVER.

%   By Andrew Kositsky
%   Copyright 2009-2010 Tectonics Observatory
%   $Revision: 1.0.0.0 $  $Date: 2009/11/30  $


%%% Verify code_dir is in fact a directory
if(~exist(code_dir,'dir'))
   error(['Directory of the Code Path ', code_dir,  ' does not exist', ...
       'Please check directory structure and spelling']); 
end

%%% Add all sub-directories of code_dir to the path recursively
add_subdir_to_path(code_dir);  % Note that code_dir must be an absolute 
                               % path to work consistently.
disp(['Subdirectories of ' code_dir ' added to path.']);