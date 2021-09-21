% Script to create the ouput files and folders names

% Marion Thomas, Harsha Bhat, last modified June 2018

%CALLS: 

%==========================================================================
%% 
disp(' ');disp('========================================');
disp('Output folders');
disp('========================================');

%% DEFINE FOLDERS AND FILES NAMES

resfolder=['res',num2str(hD),'m'];
namefoldD=strtok(datadirD(length(dir1)+1:end-1));
namefoldE=strtok(datadirE(length(dir1)+1:end-1));

%Add the parent folders
namefoldD=[resfolder,'/',nfault,namefoldD,'/'];
namefoldE=[resfolder,'/',nfault,namefoldE,'/'];
nfault2=[resfolder,'/',nfault];


%% CREATE THE FOLDERS FOR THE OUTPUT

%If no files is save
if (save_tag<1)
    
    disp('Figures are not saved, no folders are created')
    checkvarD=0;checkvarE=0;

else
    
    if exist(resfolder) == 7;  disp(['folder ',resfolder,' already exits']); ...
    else; mkdir(resfolder); disp(['folder ',resfolder]);end

    if exist(nfault2) == 7;  disp(['folder ',nfault2,' already exits']); ...
    else; mkdir(nfault2); disp(['folder ',nfault2]);end

    if exist(namefoldD) == 7,disp(['!! Careful !! folder already exits: ',namefoldD]);checkvarD=1;
    else;mkdir(namefoldD); disp(['folder ',namefoldD]); checkvarD=0;end    

    if exist(namefoldE) == 7,disp(['!! Careful !! folder already exits: ',namefoldE]);checkvarE=1;
    else;mkdir(namefoldE); disp(['folder ',namefoldE]); checkvarE=0;end    
end

