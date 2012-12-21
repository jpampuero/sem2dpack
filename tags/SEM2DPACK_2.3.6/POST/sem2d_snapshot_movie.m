% SEM2D_SNAPSHOT_MOVIE generates an animation from SEM2DPACK snapshots
%
% SYNTAX	sem2d_snapshot_movie(fname,grid)
%		sem2d_snapshot_movie(fname,grid,'Property1',Value1,'Property2',Value2,...)
%
% INPUTS	fname	code name of the snapshot field, given by the
%			initial letters of the snapshot data file name
%			(xxx_*_sem2d.dat). See also SEM2D_SNAPHOT_READ.
%			Does not work for composite fields (dmg, pla, etc)
%		grid	spectral element grid (see SEM2D_READ_SPECGRID)
%
%		Property/Value input pairs (property names are case-insensitive):
%
%		CAxis	color scale, see SEM2D_SNAPSHOT_PLOT (default = 'minmax')
%		DataDir	name of data directory (default = current directory)
%		NCycles	number of times the movie is played (default = 3)
%		TWin  	first and last snapshot in movie (default = [1 inf], plays all)
%		TPause	pause time between snapshots, in seconds (default = 0)
% 		MakeGif	export ('on') or not ('off') an animated GIF file (default = 'off')
%		GifRate	snapshots per second in animated GIF (default = 5)
%		GifFile	name of GIF file, without the .gif suffix (default = fname)
%			
function sem2d_snapshot_movie(fname,grid,varargin)

% defaults
DATADIR = [];
NCYCLES = 3;
CAXIS = [];
TWIN = [1 inf];
TPAUSE = 0;
MAKEGIF = 'off';
GIFRATE = 5; 
GIFFILE = fname;

Parse_Inputs(varargin{:});

GIFFILE = [GIFFILE '.gif'];
warning('off','MATLAB:writegif:loopcountInAppendMode');

figure(gcf);
set(gcf,'Color','w');
wrk = [];

warning('off','sem2d_snapshot_read:FileNotFound');

for n=1:NCYCLES,
  k=0;
  while 1
    k=k+1;
    if k<TWIN(1), continue, end
    if k>TWIN(2), break, end
    field = sem2d_snapshot_read(fname,k,DATADIR);
    if isstruct(field), break, end
    if isempty(field), break, end  % if last snapshot
    clf
    wrk=sem2d_snapshot_plot(field,grid,CAXIS,wrk);
    title( [fname ' snapshot # ' num2str(k)])
    drawnow
    pause(TPAUSE)

    if strcmp(lower(MAKEGIF),'on') & n==1
      I=getframe(gcf);
      I = frame2im(I);
      [X, map] = rgb2ind(I, 128);
      if k==1, wm='overwrite'; else, wm = 'append'; end
      imwrite(X,map,GIFFILE,'GIF','WriteMode',wm,'DelayTime',1/GIFRATE,'LoopCount',NCYCLES);  
    end

  end
  if k==TWIN(1), break, end % if field does not exist
end

warning('on','MATLAB:writegif:loopcountInAppendMode');
warning('on','sem2d_snapshot_read:FileNotFound');

%-----------
% PARSE_INPUTS sets variables in the caller function according to inputs.
% 	Override defaults by Property/Value pairs
%       By convention the variables to set have UPPER CASE names in the caller. 
function Parse_Inputs(varargin)

if isempty(varargin), return, end
for k=2:2:length(varargin),
  assignin('caller', upper(varargin{k-1}), varargin{k} );
end
