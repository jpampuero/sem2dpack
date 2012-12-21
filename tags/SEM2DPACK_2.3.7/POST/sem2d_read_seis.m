% SEM2D_READ_SEIS reads seismogram outputs from SEM2DPACK
%
% SYNTAX	data = sem2d_read_seis()
% 		data = sem2d_read_seis(dirname)
%
% INPUTS	dirname	[pwd]	data directory
%
% OUTPUTS	data	structure containing:
%			nt	number of timesteps
%			dt	timestep
%			nsta	number of stations
%			x,z	position of stations
%			ux,uy,uz	seismograms, size=[nt,nsta]
%				
function data = sem2d_read_seis(dirname)

if exist('dirname','var')
  dir1=pwd;
  cd(dirname);
end

% Read parameters from header file
[data.dt,data.nt,data.nsta] = textread('SeisHeader_sem2d.hdr','%n%n%n',1,'headerlines',1);
[data.x,data.z] = textread('SeisHeader_sem2d.hdr','%f%f','headerlines',3);

% Read seismograms
usize = [data.nt,data.nsta];
fid=fopen('Ux_sem2d.dat'); 
if fid ~= -1
  data.ux = fread(fid,usize,'single') ; 
  fclose(fid);
end
fid=fopen('Uy_sem2d.dat'); 
if fid ~= -1
  data.uy = fread(fid,usize,'single') ; 
  fclose(fid);
end
fid=fopen('Uz_sem2d.dat'); 
if fid ~= -1
  data.uz = fread(fid,usize,'single') ; 
  fclose(fid);
end

if exist('dirname','var'), cd(dir1); end
