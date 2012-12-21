% SEM2D_SNAPSHOT_READ reads a snapshot output from SEM2DPACK
%
% SYNTAX	field = sem2d_snapshot_read(fname,isnap,datadir)
%
% INPUTS	fname	code name of the snapshot field, given by the
%			initial letters of the snapshot data file name
%			(xxx_*_sem2d.dat):
%
%   	                Field \ Type	P-SV			SH
%			------------------------------------------------
%       		displacement    dx,dz                   dy
%       		velocity        vx,vz                   vy
%       		acceleration    ax,az                   ay
%       		strain          e11,e22,e12             e13,e23
%       		stress          s11,s22,s12 	        s13,s23
%			divergence	div			--
%			curl		curl			--
%			damage		dmg			--
%			plastic		pla			--
%
% 		isnap   snapshot index, given by the numbers of the 
%			snapshot data file name (*_xxx_sem2d.dat)
%		datadir	data directory name (default: current directory)
%
% OUTPUT	field	array (or structure of arrays) containing snapshot data, 
%			in either of two formats:
% 			 1. global format, node-by-node (d, v, a)
% 			    ex: if fname='vx', field(k) is the x-component of velocity
%			    of the k-th global node
% 			 2. local format, element-by-element (e, s, dmg)
% 			    ex: if fname='e11', field(i,j,e) is the xx component of
%			    strain at the (i,j)-th local node of the e-th element
%			If the snapshot data file contains multiple variables, 
%			field is a structure of arrays. For instance, if fname='dmg'
%			field contains the damage variable (alpha) and the plastic
%			strains (ep11, ep22, ep12)
%
function field = sem2d_snapshot_read(fname,isnap,datadir)

if exist('isnap','var') && ~isempty(isnap)
  fname = sprintf('%s_%3.3u_sem2d.dat',fname,isnap);
else
  fname = sprintf('%s_sem2d.dat',fname);
end
if exist('datadir','var') && ~isempty(datadir)
  if datadir(end) ~= '/', datadir = [datadir '/']; end
else
  datadir='';
end
fname = [datadir fname];

if ~exist(fname,'file'),
  field = [];
  warning('sem2d_snapshot_read:FileNotFound',['File ' fname ' not found'])
  return 
end

fid = fopen(fname); 
field = fread(fid,inf,'single'); 
fclose(fid);

[nelem,npoin,ngll] = textread([datadir 'grid_sem2d.hdr'],'%u%*u%*u%u%u','headerlines',1);

% if reading damage,plastic,etc variables 
% create a structure and reshape fields in element-wise storage
fname_pre = filename_prefix(fname);
if ~isempty(strfind(fname_pre,'dmg_'))
 % WARNING : should read dmg_elems.tab 
  ftmp = reshape( field, ngll,ngll,[],nelem );
  field = [];
  field.alpha = squeeze(ftmp(:,:,1,:)); % damage variable
  field.ep11 = squeeze(ftmp(:,:,2,:));  % plastic strain
  field.ep22 = squeeze(ftmp(:,:,3,:));
  field.ep12 = squeeze(ftmp(:,:,4,:));

elseif ~isempty(strfind(fname_pre,'pla_'))
 % WARNING : should read pla_elems.tab 
  ftmp = reshape( field, ngll,ngll,[],nelem );
  field = [];
  field.ep11 = squeeze(ftmp(:,:,1,:));  % plastic strain
  field.ep22 = squeeze(ftmp(:,:,2,:));
  field.ep12 = squeeze(ftmp(:,:,3,:));

% reshape fields in element-wise storage
elseif length(field)~=npoin 
  field=reshape(field,[ngll,ngll,nelem]); 
end


%----
% extracts the prefix from a file name 
% ex: filename_prefix('/home/me/toto/stuff.dat') gives 'stuff'
%
function pre = filename_prefix(fname)

k1= strfind(fname,'/');
if isempty(k1), k1=0; end
k1=k1(end)+1;

k2= strfind(fname,'.');
if isempty(k2), k2=length(fname)+1; end 
k2=k2(end)-1;

pre = fname(k1:k2);
