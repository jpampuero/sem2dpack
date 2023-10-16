% SEM2D_READ_FAULT reads fault outputs from SEM2DPACK
%
% SYNTAX	data = sem2d_read_fault(name)
%
% INPUT		name	[Flt*] 	prefix of header and data files (name_sem2d.*) 
%				The default is the first FltXX_sem2d.* found 
%				in the current directory.
%
% OUTPUT	data.nx		number of fault nodes
%		data.nt		number of time samples
%		data.dt		time step
%		data.x,data.z	coordinates of fault nodes
%		data.d		slip [nx,nt]
%		data.v		slip rate
%		data.st		shear stress
%		data.sn		normal stress
%		data.mu		friction coefficient
%		data.theta  state variable	
%		data.st0	initial value of shear stress
%		data.sn0	initial value of normal stress
%		data.mu0	initial value of friction coefficient
%		data.theta0	initial value of state variable 
% 
%		If output on each side of the fault (osides=T):
%  		data.d1t	displacement on side 1, fault parallel component
%  		data.d2t	displacement on side 2, fault parallel component
%  		data.d1n	displacement on side 1, fault normal component
%  		data.d2n	displacement on side 2, fault normal component
%  		data.v1t	velocity on side 1, fault parallel component
%  		data.v2t	velocity on side 2, fault parallel component
%  		data.v1n	velocity on side 1, fault normal component
%  		data.v2n	velocity on side 2, fault normal component
%
% NOTE		Fault normal components on each side of the fault are exported only in P-SV
%
function data = sem2d_read_fault(name, wdir)

% length of the tag at the begining and end of a binary record
% in number of single precision words (4*bytes)
LENTAG = 2; % gfortran older versions
LENTAG = 1;

if nargin<2
    wdir=pwd;
end

% assumes header file name is FltXX_sem2d.hdr
if ~exist('name','var')
  list = dir('Flt*.hdr');
  list = {list.name};
  if isempty(list)
    name = '';
  else
    name=list{1}(1:5);
  end
end

% Read parameters from header file
hdr = strcat(wdir,'/',name,'_sem2d.hdr');
if ~exist(hdr,'file')
  data=[]; 
  warning(['File ' hdr ' not found'])
  return
end

% check if the time info file exist or not
tinfo = strcat(wdir,'/',name,'_time_sem2d.tab');
if exist(tinfo, 'file')
    %adaptive time stepping from rate state cycle simulation
    fprintf('Variable time step...\n');
    [data.it, data.dt, data.t, data.eqnum, data.isdynamic, ...
    data.isswitch, data.isEQ]=sem2d_read_rsf_timeinfo(tinfo);
    data.nt = length(data.t);
    [data.nx,ndat] = textread(hdr,'%n%n',1,'headerlines',1);
else
    fprintf('Fixed time step...\n');
    [data.nx,ndat,data.nt,data.dt] = textread(hdr,'%n%n%n%n',1,'headerlines',1);
end

[data.x,data.z] = textread(hdr,'%f%f','headerlines',4);

% Read initial fault data
if exist([wdir,'/',name '_init_sem2d.tab'],'file')
  raw = load([wdir,'/',name '_init_sem2d.tab']);
  data.st0 = raw(:,1);
  data.sn0 = raw(:,2);
  data.mu0 = raw(:,3);
  if size(raw, 2)>3
     data.phi   = raw(:, 4);% Marco's implemention added phi
     data.theta0=raw(:,5); 
     data.v0=raw(:,6); 
     data.rsf_a = raw(:, 7);
     data.rsf_b = raw(:, 8);
  end
end

% Read fault data in a big matrix
dat  = strcat(wdir, '/', name,'_sem2d.dat');
fid=fopen(dat); 
raw = fread(fid,[data.nx+2*LENTAG,inf],'single') ; 
fclose(fid);
%raw = reshape(raw(2:data.nx+1,:),[data.nx ndat data.nt]);
raw = reshape(raw(LENTAG+1:LENTAG+data.nx,:),data.nx,ndat,[]);

% Reformat each field [nx,nt]
data.d  = squeeze(raw(:,1,:)); 
data.v  = squeeze(raw(:,2,:)); 
data.st = squeeze(raw(:,3,:)); 
data.sn = squeeze(raw(:,4,:)); 
data.mu  = squeeze(raw(:,5,:)); 
switch ndat
    case 5+1
      data.theta=squeeze(raw(:,6,:));
    case 5+2
      data.theta=squeeze(raw(:,6,:));
      data.phi = squeeze(raw(:,7,:));
    case 5+4
      data.d1t  = squeeze(raw(:,6,:)); 
      data.d2t  = squeeze(raw(:,7,:)); 
      data.v1t  = squeeze(raw(:,8,:)); 
      data.v2t  = squeeze(raw(:,9,:)); 
    case 5+4*2
      data.d1t  = squeeze(raw(:,6,:)); 
      data.d1n  = squeeze(raw(:,7,:)); 
      data.d2t  = squeeze(raw(:,8,:)); 
      data.d2n  = squeeze(raw(:,9,:)); 
      data.v1t  = squeeze(raw(:,10,:)); 
      data.v1n  = squeeze(raw(:,11,:)); 
      data.v2t  = squeeze(raw(:,12,:)); 
      data.v2n  = squeeze(raw(:,13,:)); 
    case 5+4+1
      data.d1t  = squeeze(raw(:,6,:)); 
      data.d2t  = squeeze(raw(:,7,:)); 
      data.v1t  = squeeze(raw(:,8,:)); 
      data.v2t  = squeeze(raw(:,9,:)); 
      data.theta  = squeeze(raw(:,10,:)); 
    case 5+4*2+1
      data.d1t  = squeeze(raw(:,6,:)); 
      data.d1n  = squeeze(raw(:,7,:)); 
      data.d2t  = squeeze(raw(:,8,:)); 
      data.d2n  = squeeze(raw(:,9,:)); 
      data.v1t  = squeeze(raw(:,10,:)); 
      data.v1n  = squeeze(raw(:,11,:)); 
      data.v2t  = squeeze(raw(:,12,:)); 
      data.v2n  = squeeze(raw(:,13,:)); 
      data.theta  = squeeze(raw(:,14,:)); 
end
end

