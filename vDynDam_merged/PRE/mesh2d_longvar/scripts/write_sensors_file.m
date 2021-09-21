% Script to create the ouput file with the sensors location

% Marion Thomas, last modified November 2018

%CALLS: 

%==========================================================================
%% CORRECTIONS BASED ON FAULT COORDINATES

%find the unique x coordinates 
xloc1=unique(xsens);

%only keep the ones that correspond to a x coordinates one the fault
xloc=xloc1(find(xloc1>=xF(1) & xloc1 <= xF(end)));

%find the corresponding z coordinates on the fault
z_xloc=zeros(size(xloc));
for i=1:numel(xloc)
    dxloc=abs(xloc(i)-xF);
    z_xloc(i)=zF(find(dxloc==min(dxloc)));
end

%Define the z coordinates of the sismometers
zsens=dzsens;
for i=1:numel(xloc)
    idS=find(xsens == xloc(i));
    zsens(idS(1:end-1))=zsens(idS(1:end-1))+z_xloc(i);
end

%% PLOT TO CHECK POSITIONS

figure(2); hold on
plot(xsens,zsens,'r*')

%% WRITE SENSORS FILE

%file name
sensorfile='sensors.ini';

%Damage folder
wfile = [rootD,namefoldD,'/',sensorfile];
fid = fopen(wfile,'wt');
fprintf(fid,'%6.2f %12.0f\r\n',[xsens zsens]');
fprintf(fid,'%6.0f %12.0f\r\n',[nuclocX zF(xF==nuclocX)]');
fclose(fid);

%Elastic folder
if nodmg_tag == 1
wfile2 = [rootD,namefoldE,'/',sensorfile];
fid2 = fopen(wfile2,'wt');
fprintf(fid2,'%6.0f %12.0f\r\n',[xsens zsens]');
fprintf(fid2,'%6.0f %12.0f\r\n',[nuclocX zF(xF==nuclocX)]');
fclose(fid2);
end