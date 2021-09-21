% ************************************************************************
% This function reads the info about simulations from the file
% "sem2dpack.sta"

% Marion Thomas, last modified July 2018
% CALLS: 

%% ************************************************************************
function [info,MuS,MuD,ndof,Nmat] = getsimulinfo(datadir)

fid = fopen([datadir 'sem2dpack.sta'], 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

%% Numbers of materials

line = find(ismember(C{1}(:),'M a t e r i a l   P r o p e r t i e s'));
str = char(C{1}(line+3)); 
disp(strrep(strrep(strrep(str, '. ', '.'),'.',''),'(',' ('))
% disp(' ')
ind = strfind(str,'=');
Nmat=str2double(strtrim(str(ind+1:end)));
info=zeros(Nmat,5);

%% FIND THE BULK PROPERTIES

for i=1:Nmat
    
    %find the right material
    tag=num2str(i);
    disp(['#---- Material ',tag,' -----------------------'])
    line = find(ismember(C{1}(:),['Material index. . . . . . . . . . . (tag) = ',tag]));

    %find value of initial damage
    str = char(C{1}(line+21)); 
    ind = strfind(str,'=');
    textdisp=strrep(strrep(strrep(strrep(str(1:ind-1), '. ', '.'),'Initial ',''),'.',''),'(',' (');
    info(i,1) = str2double(strtrim(str(ind+1:end)));
    disp([textdisp,'= ',num2str(info(i,1))])

    %find value of initial cs
    str = char(C{1}(line+12)); 
    ind = strfind(str,'=');
    textdisp=strrep(strrep(strrep(str(1:ind-1), '. ', '.'),'.',''),'(',' (');
    info(i,2) = str2double(strtrim(str(ind+1:end)));
    disp([textdisp,'= ',num2str(info(i,2)),' m/s'])

    %find value of initial cp
    str = char(C{1}(line+11)); 
    ind = strfind(str,'=');
    textdisp=strrep(strrep(strrep(str(1:ind-1), '. ', '.'),'.',''),'(',' (');
    info(i,3) = str2double(strtrim(str(ind+1:end)));
    disp([textdisp,'= ',num2str(info(i,3)),' m/s'])
    
    %find value of initial density rho
    str = char(C{1}(line+14)); 
    ind = strfind(str,'=');
    textdisp=strrep(strrep(strrep(str(1:ind-1), '. ', '.'),'.',''),'(',' (');
    info(i,4) = str2double(strtrim(str(ind+1:end)));
    disp([textdisp,'   = ',num2str(info(i,4))])

    %find value of initial NU
    str = char(C{1}(line+15));
    ind = strfind(str,'=');
    textdisp=strrep(strrep(strrep(str(1:ind-1), '. ', '.'),'.',''),'(',' (');
    info(i,5)= str2double(strtrim(str(ind+1:end)));
    disp([textdisp,'= ',num2str(info(i,5))])
    
end

%% FIND THE FRICTION PROPERTIES

disp(' '); 
disp('========================================');
disp('Fault Parameters')
disp('========================================');

%find value of initial Mus
line = find(ismember(C{1}(:),'B o u n d a r y   C o n d i t i o n s'));
str = char(C{1}(line+25)); 
disp(strrep(strrep(strrep(str, '. ', '.'),'.',''),'(',' ('))
ind = strfind(str,'=');
MuS = str2double(strtrim(str(ind+1:end)));

%find value of initial MuD
str = char(C{1}(line+26)); 
disp(strrep(strrep(strrep(str, '. ', '.'),'.',''),'(',' ('))
ind = strfind(str,'=');
MuD = str2double(strtrim(str(ind+1:end)));

%find value of ndof
line = find(ismember(C{1}(:),'G e n e r a l   P a r a m e t e r s'));
str = char(C{1}(line+5));
disp(strrep(strrep(strrep(str, '. ', '.'),'.',''),'(',' ('))
ind = strfind(str,'=');
ndof = str2double(strtrim(str(ind+1:end)));


