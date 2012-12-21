% MESH2D_READ 	reads a 2D mesh database from a *.mesh2d file 
%
% SYNTAX	mesh = mesh2d_read(fname)
% 
% INPUTS	fname	mesh file name 
%			The suffix .mesh2d is automatically added to the file name, if absent.
%
% OUTPUTS	mesh	a mesh2d structure (see MESH2D_TFI)
%
function m=mesh2d_read(fname)

fname = deblank(fname);
if isempty(regexp(fname,'\.mesh2d$')), fname = strcat(fname,'.mesh2d'); end

fid = fopen(fname,'r');
line = fgetl(fid);
dat=fscanf(fid,'%u %u %u %u\n');
nel=dat(1);
q=dat(2);
np=dat(3);
nb=dat(4);
line = fgetl(fid);
m.coor = fscanf(fid,'%*u %f %f\n',[2 np]);
line = fgetl(fid);
fmt = ['%*u ' repmat('%u ',1,q+1) '\n'];
dat = fscanf(fid,fmt,[q+1 nel]);
m.enod = dat(1:q,:);
m.etag = dat(q+1,:);
for k=1:nb,
  line = fgetl(fid);
  dat=fscanf(fid,'%u %u\n');
  bctag = dat(1);
  nbel = dat(2);
  line = fgetl(fid);
  m.bnds{bctag} = fscanf(fid,'%*u %u %u\n',[2 nbel]);
end
fclose(fid);

