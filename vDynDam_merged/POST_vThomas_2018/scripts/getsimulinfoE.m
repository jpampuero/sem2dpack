% ************************************************************************
% This function reads the info about simulations from the file
% "sem2dpack.sta"

% Marion Thomas, last modified July 2018
% CALLS: 

%% ************************************************************************
function [info,Nmat] = getsimulinfo(datadir)

fid = fopen([datadir 'sem2dpack.sta'], 'r');
C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

%% Numbers of materials

line = find(ismember(C{1}(:),'M a t e r i a l   P r o p e r t i e s'));
str = char(C{1}(line+3)); 
disp(strrep(strrep(strrep(str, '. ', '.'),'.',''),'(',' ('))
ind = strfind(str,'=');
Nmat=str2double(strtrim(str(ind+1:end)));
info=zeros(Nmat,3);

%% FIND THE BULK PROPERTIES

for i=1:Nmat
    
    %find the right material
    tag=num2str(i);
    disp(['#---- Material ',tag,' -----------------------'])
    line = find(ismember(C{1}(:),['Material index. . . . . . . . . . . (tag) = ',tag]));

   %find value of initial cs
    str = char(C{1}(line+3)); 
    disp(strrep(strrep(strrep(str, '. ', '.'),'.',''),'(',' ('))
    ind = strfind(str,'=');
    info(i,2) = str2double(strtrim(str(ind+1:end)));

    %find value of initial cp
    str = char(C{1}(line+2)); 
    disp(strrep(strrep(strrep(str, '. ', '.'),'.',''),'(',' ('))
    ind = strfind(str,'=');
    info(i,3) = str2double(strtrim(str(ind+1:end)));

    %find value of initial density
    str = char(C{1}(line+4)); 
    disp(strrep(strrep(strrep(str, '. ', '.'),'.',''),'(',' ('))
    ind = strfind(str,'=');
    info(i,4) = str2double(strtrim(str(ind+1:end)));
end


