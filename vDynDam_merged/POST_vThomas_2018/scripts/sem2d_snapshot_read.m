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
if ~isempty(strfind(fname_pre,'elm_'))
 % WARNING : should read pla_elems.tab 
  ftmp = reshape( field, ngll,ngll,[],nelem );
  sf=size(ftmp,3);
  field = [];
  field.alpha = squeeze(ftmp(:,:,1,:)); % state variable
  if sf == 3
      field.invI = squeeze(ftmp(:,:,2,:));  % strain
      field.invII = squeeze(ftmp(:,:,3,:));
  end
  
elseif ~isempty(strfind(fname_pre,'dmg_'))
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

%Dynamic Damage
elseif ~isempty(strfind(fname_pre,'dyD_'))
  ftmp = reshape( field, ngll,ngll,[],nelem );
  sf=size(ftmp,3);
  field = [];
  field.D = squeeze(ftmp(:,:,1,:)); % state variable
  if sf >= 3
      field.cs = squeeze(ftmp(:,:,2,:));  % Apparent S and P waves
      field.cp = squeeze(ftmp(:,:,3,:));
  end
  if sf >= 4
      field.dldt = squeeze(ftmp(:,:,4,:));% wing-crack speed
  end
  if sf >= 5
      field.R = squeeze(ftmp(:,:,5,:));  % Regime
  end
  if sf >= 7
      field.invI = squeeze(ftmp(:,:,6,:));  % Apparent S and P waves
      field.invII = squeeze(ftmp(:,:,7,:));
  end
  if sf >= 8
      field.Rcum = squeeze(ftmp(:,:,8,:));  % stress intensity factor
  end
  if sf >= 11
      field.A = squeeze(ftmp(:,:,9,:));  % KI parameters
      field.B = squeeze(ftmp(:,:,10,:));
      field.C = squeeze(ftmp(:,:,11,:));
  end

elseif length(field)~=npoin 
  field=reshape(field,[ngll,ngll,nelem]); 
end
%   integer, intent(in) :: ndof
%   type(matwrk_dyndmg_type), intent(in) :: m
%   real :: dat(size(m%D,1),size(m%D,2),4)
%   dat(:,:,1) = real(m%D)
%   dat(:,:,2) = real(m%css)
%   dat(:,:,3) = real(m%cps)
%   dat(:,:,4) = real(m%v)
% !  dat(:,:,5) = real(m%KI)
% !  dat(:,:,6) = real(m%A)
% !  dat(:,:,7) = real(m%B)
% !  dat(:,:,8) = real(m%C)
% !  dat(:,:,9) = real(m%R)







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
