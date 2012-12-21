% sem2d_view(scale)
%
% PURPOSE	interactively plots color snapshots 
%		of data output from SEM2DPACK
%		
% INPUT		scale	[min max] color scale bounds 
%			The default [] sets automatic scale
%
function sem2d_view(scale)

if ~exist('scale','var'), scale=[]; end

% list available output data files
% data file name format: xx_###_sem2d.dat
s=unix('ls -1 *_[0-9]*_sem2d.dat > tmp');
if s, return, end 
list=textread('tmp','%s');

% cut the field names
showlist=cell2mat(list);
showlist=showlist(:,1:7);

%-- Prepare grid data :

% Grid parameters :
[nelem,npgeo,ngnod,npoin,ngll] = ...
  textread('grid_sem2d.hdr','%u%u%u%u%u','headerlines',1);

% coord(npoin,2) = coordinates of the GLL nodes
fid=fopen('coord_sem2d.dat');
coord = fread(fid,[2,npoin],'single')' ; 
fclose(fid);
fid=fopen('ibool_sem2d.dat'); 

% ibool(ngll,ngll,nelem) = local to global numbering
ibool = fread(fid,inf,'int'); 
fclose(fid);
ibool = reshape(ibool,[ngll,ngll,nelem]);

% cell data
indx_nodal = Init2dSnapshot(ibool);
[indx_elem,coord_elem] = Init2dSnapshot(ibool,coord);

%-- Interactive plots:

% Window 1:	menu 1 --> select field and snapshot
%	   	menu 2 --> select scale
% Window 2: 	plot 
k=1;
while 1

  [k,ok] = listdlg('ListString',showlist,'InitialValue',k,...
                   'SelectionMode','single', ...
                   'Name','SEM2DPACK view fields' , ...
                   'PromptString','FieldCode_SnapshotIndex');
  if ~ok, break, end

  fname = list{k};
  fid = fopen(fname); field = fread(fid,inf,'single'); fclose(fid);
  clf
  if isempty(findstr(fname(1),'es'))  % is nodal
    Plot2dSnapshot(field,coord,indx_nodal,scale);
  else
    %field=reshape(field,[ngll,ngll,nelem]);
    Plot2dSnapshot(field,coord_elem,indx_elem,scale);
  end

end
