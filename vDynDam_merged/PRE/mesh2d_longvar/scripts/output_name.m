% Script to create the ouput files and folders names

% Marion Thomas, Harsha Bhat, last modified June 2018

%CALLS: 

%==========================================================================
%% 

disp(' ');disp('========================================');
disp(''); disp('Output folders');disp(' ')

%% DEFINE FOLDERS AND FILES NAMES

%resolution
resfolder=['res',num2str(h),'m'];

%Type of fault
if exist('rms')==1
    namefoldD=['dydmg_roughF_' ,num2str(s.Seed)];
    namefoldE=['nodmg_roughF_' ,num2str(s.Seed)];
    namemesh=[nfault,'_',num2str(s.Seed)];
else
    namefoldD='dydmg';
    namefoldE='nodmg';
    namemesh=nfault;
end

%D0 distribution
if own_tag==1
    namefoldD=[namefoldD,'_ownD0'];
elseif exp_tag==1
    namefoldD=[namefoldD,'_expD0'];
elseif gaus_tag==1
    namefoldD=[namefoldD,'_gausD0'];
else
    NumD0=numel(unique(D0));
    namefoldD=[namefoldD,'_',num2str(NumD0),'D0'];
end

%Bimaterial
if Nmat==2, namefoldD=['bimat_',namefoldD];namefoldE=['bimat_',namefoldE]; end

%Add the parent folders
namefoldD=[resfolder,'/',nfault,'/',namefoldD];
namefoldE=[resfolder,'/',nfault,'/',namefoldE];
nfault2=[resfolder,'/',nfault];

%% CREATE THE FOLDERS FOR THE OUTPUT

if exist(resfolder) == 7;  disp(['folder ',resfolder,' already exits']); ...
else; mkdir(resfolder); disp(['folder ',resfolder]);end

if exist(nfault2) == 7;  disp(['folder ',nfault2,' already exits']); ...
else; mkdir(nfault2); disp(['folder ',nfault2]);end

if exist(namefoldD) == 7,disp(['!! Careful !! folder already exits: ',namefoldD]);checkvarD=1;
else;mkdir(namefoldD); disp(['folder ',namefoldD]); checkvarD=0;end

if nodmg_tag == 1
if exist(namefoldE) == 7, disp(['!! Careful !! folder already exits: ',namefoldE]);checkvarE=1; ...
else; mkdir(namefoldE); disp(['folder ',namefoldE]);checkvarE=0;end
end



%% But need to be on the main script

%     while(1)
%         m=input('Do you want to continue, Y/N [Y]:','s');
%         if m=='N'; break; end
%         if m=='n'; break; end
%     end
% while(1)
% 
% while checkvar==1
%     promptMessage = sprintf('Do you want to Continue processing,\nor Cancel to abort processing?');
%     button = questdlg(promptMessage, 'Continue', 'Continue', 'Cancel', 'Continue');
%     if strcmpi(button, 'Cancel')
%         return; % Or break or continue
%     end
% %         m=input('Do you want to continue, Y/N [Y]:','s');
% %         if strcmp(m,'N')==1;'coucou'; return; end
% %         if strcmp(m,'n')==1;'coucou'; return; end
% 
% end
