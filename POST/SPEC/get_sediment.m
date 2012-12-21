% This script gets the sediment thickness under each station
% INPUT: datapath,sedimentfile
% OUTPUT: thickness
% sedimentfile must contain ALL basement points
% including the closing points at the edges of the valley !

var1 = load (sedimentfile) ; 
cd(datapath); var2 = load('xsismos.dat');cd (sourcepath);

size(var1);
N1 = ans(1);
size(var2);
N2 = ans(1);

col_x = 1;
col_z = 2;

new = interp1( var1(1:N1,col_x), var1(1:N1,col_z), var2(1:N2,col_x), 'linear');
var2(1:N2,col_z) = var2(1:N2,col_z) - new;

% set to zero depth the stations on rock:
for i=1:N2
    if isnan(var2(i,col_z))
        var2(i,col_z) = 0.;
    end
end

%plot(var1(1:N1,col_x), var1(1:N1,col_z),   var2(1:N2,col_x),var2(1:N2,col_z));
%save -ascii Depthfile var2

thickness = var2 ;
