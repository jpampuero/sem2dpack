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
%   data.sts	shear stress T "stick" (for rupture velocity calculations)
%		data.st0	initial value of shear stress
%		data.sn0	initial value of normal stress
%		data.mu0	initial value of friction coefficient
%   data.wei	node weights
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
function data = sem2d_read_fault(name)

% length of the tag at the begining and end of a binary record
% in number of single precision words (4*bytes)
LENTAG = 2; % gfortran older versions
LENTAG = 1;

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
hdr = strcat(name,'_sem2d.hdr');
if ~exist(hdr,'file')
  data=[]; 
  warning(['File ' hdr ' not found'])
  return
end
[data.nx,ndat,data.nt,data.dt] = textread(hdr,'%n%n%n%n',1,'headerlines',1);
[data.x,data.z] = textread(hdr,'%f%f','headerlines',4);

% Read initial fault data
if exist([name '_init_sem2d.tab'],'file')
  raw = load([name '_init_sem2d.tab']);
  data.st0 = raw(:,1);
  data.sn0 = raw(:,2);
  data.mu0 = raw(:,3);
  if size(raw,2) == 4
      data.wei = raw(:,4);            % save weights if they are in the file
  end
end

% Read fault data in a big matrix
dat  = strcat(name,'_sem2d.dat');
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

% Backwards compatibility with output before Tstick
if ismember(ndat, [6, 6+4, 6+4*2])
    data.sts = squeeze(raw(:,6,:)); 
    core = 6;
elseif ismember(ndat, [5, 5+4, 5+4*2])
    core = 5;
end

if ismember(ndat, [6+4, 5+4])
  data.d1t  = squeeze(raw(:,core+1,:)); 
  data.d2t  = squeeze(raw(:,core+2,:)); 
  data.v1t  = squeeze(raw(:,core+3,:)); 
  data.v2t  = squeeze(raw(:,core+4,:)); 
elseif ismember(ndat, [6+4*2, 5+4*2])
  data.d1t  = squeeze(raw(:,core+1,:)); 
  data.d1n  = squeeze(raw(:,core+2,:)); 
  data.d2t  = squeeze(raw(:,core+3,:)); 
  data.d2n  = squeeze(raw(:,core+4,:)); 
  data.v1t  = squeeze(raw(:,core+5,:)); 
  data.v1n  = squeeze(raw(:,core+6,:)); 
  data.v2t  = squeeze(raw(:,core+7,:)); 
  data.v2n  = squeeze(raw(:,core+8,:)); 
end
