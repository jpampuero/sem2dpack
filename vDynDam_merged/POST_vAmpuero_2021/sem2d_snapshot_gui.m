% SEM2D_SNAPSHOT_GUI interactive plot of snapshot outputs from SEM2DPACK
%
% SYNTAX	sem2d_snapshot_gui(datadir)
%		
% INPUT		datadir	name of the directory containing the SEM2DPACK output
%		        (default is the current directory)
%
function sem2d_snapshot_gui(inputdir)

if exist('inputdir','var'), 
  if ~isdir(inputdir)
    my_error('Invalid directory name'); 
    defdir = pwd;
  else
    defdir=fix_path(inputdir);
  end
else
  defdir = pwd;
end

% defaults
data.dir = defdir;
data.times = [];
data.fields = {''};
data.subfields = {''};
data.scale = 'minmax';

f = figure('Visible','off','Position',[0 300 300 400], ...
           'Name','SEM2DPACK snapshot viewer','NumberTitle','off', ...
           'MenuBar','none','ToolBar','none','Units','normalized');

% Select directory
dir_panel = uipanel('Title','Data directory','Position',[0.02 0.85 0.96 0.13]);
dir_edit = uicontrol(dir_panel,'Style','edit', ...
                     'Units','normalized','Position',[0.02 0.02 0.58 0.96], ...
                     'String',data.dir,...
		     'Callback',{@dir_edit_callback});
dir_browse = uicontrol(dir_panel,'Style','pushbutton', ...
                       'Units','normalized','Position',[0.62 0.02 0.17 0.96],...
                       'String','Browse',...
		       'Callback',{@dir_browse_callback});
dir_load = uicontrol(dir_panel,'Style','pushbutton', ...
                     'Units','normalized','Position',[0.81 0.02 0.17 0.96],...
                     'String','Load', ...
		     'Callback',{@dir_load_callback});

% Select field
data_panel = uipanel('Title','Data','Position',[0.02 0.3 0.96 0.55]);
data_list1 = uicontrol(data_panel,'Style','listbox', ...
                       'Units','normalized','Position',[0.02 0.02 0.3 0.90], ...
		       'Value',2, ...
		       'Callback',{@data_list1_callback});
data_title1 = uicontrol(data_panel,'Style','text', ...
                        'Units','normalized','Position',[0.02 0.93 0.3 0.05],...
			'String','Time index'); 
data_list2 = uicontrol(data_panel,'Style','listbox', ...
                       'Units','normalized','Position',[0.34 0.02 0.31 0.90], ...
		       'Callback',{@data_list2_callback});
data_title2 = uicontrol(data_panel,'Style','text', ...
                        'Units','normalized','Position',[0.34 0.93 0.31 0.05],...
			'String','Field'); 
data_list3 = uicontrol(data_panel,'Style','listbox', ...
                       'Units','normalized','Position',[0.67 0.02 0.31 0.90], ...
		       'Callback',{@data_list3_callback});
data_title3 = uicontrol(data_panel,'Style','text', ...
                        'Units','normalized','Position',[0.67 0.93 0.31 0.05],...
			'String','Sub-field'); 

% Select color scale
color_panel = uipanel('Title','Color scale','Position',[0.02 0.16 0.96 0.14]);
color_menu = uicontrol(color_panel,'Style','popupmenu', ...
                       'Units','normalized','Position',[0.02 0.02 0.3 0.96], ...
		       'String',{'min - max','-/+ max abs','manual'}, ...
		       'Callback',{@color_menu_callback});
color_range = uicontrol(color_panel,'Style','text', ...
                      'Units','normalized','Position',[0.34 0.02 0.16 0.96], ...
		      'HorizontalAlignment','right', ...
		      'String','Range');
color_min = uicontrol(color_panel,'Style','edit', ...
                      'Units','normalized','Position',[0.52 0.02 0.22 0.96],'String','-1', ...
		      'Enable','off', ...
		      'TooltipString','Low bound of color scale range', ...
		      'Callback',{@color_min_callback});
color_max = uicontrol(color_panel,'Style','edit', ...
                      'Units','normalized','Position',[0.76 0.02 0.22 0.96],'String','1', ...
		      'Enable','off', ...
		      'TooltipString','High bound of color scale range', ...
		      'Callback',{@color_max_callback});

% Plot
plot_button = uicontrol(f,'Style','pushbutton', ...
                        'Units','normalized','Position',[0.02 0.02 0.47 0.12], ...
			'String','Plot', ...
		        'Callback',{@plot_callback});
close_button = uicontrol(f,'Style','pushbutton', ...
                        'Units','normalized','Position',[0.51 0.02 0.47 0.12], ...
			'String','Close', ...
		        'Callback',{@close_callback});

data_refresh;

% Move the GUI to the center of the screen.
movegui(f,'northwest')
% Make the GUI visible on screen.
set(f,'Visible','on');
% Make the GUI invisible to command line manipulations.
set(f,'HandleVisibility','callback');

%-------------------------------
function dir_edit_callback(source,eventdata)
data.dir = get(source, 'String');
data_refresh;
end

%-------------------------------
function dir_browse_callback(source,eventdata)
data.dir = uigetdir(defdir);
% if user presses 'Cancel' data_dir = 0
if isstr(data.dir)  
  set(dir_edit,'String',data.dir);
  data_refresh;
else
  data.dir = defdir;
end
end

%-------------------------------
function dir_load_callback(source,eventdata)
data_refresh;
end

%-------------------------------
function data_list1_callback(source,eventdata)
data.time = data.times(get(source,'Value'));
data_subfields_refresh;
end

%-------------------------------
function data_list2_callback(source,eventdata)
data.field = data.fields{get(source,'Value')};
data_subfields_refresh;
end

%-------------------------------
function data_list3_callback(source,eventdata)
data.subfield = data.subfields{get(source,'Value')};
end

%-------------------------------
function data_refresh

if ~isdir(data.dir)
  my_error('Invalid directory name')
  return
end

% list all data files
d=dir([data.dir '/*_*_sem2d.dat']);

if isempty(d)
  my_error('Directory does not contain SEM2D data')
  return
end

% load grid data
data.grid = sem2d_read_specgrid(data.dir);

% extract times and field names
fields = cell(length(d),1);
times = zeros(length(d),1);
for k = 1:length(d),
  loc = findstr('_',d(k).name);
  fields{k} = d(k).name(1:loc(1)-1);
  time_chk = str2num( d(k).name(loc(1)+1:loc(2)-1) );
  if isempty(time_chk) 
    times(k) = nan; 
  else
    times(k) = time_chk; 
  end
end
data.times = unique(times(isfinite(times)));
data.fields = unique(fields);

set(data_list1,'String',cellstr(num2str(data.times)));
set(data_list1,'Value', min( length(get(data_list1,'String')), get(data_list1,'Value') ) )
data.time = data.times(get(data_list1,'Value'));

set(data_list2,'String',data.fields);
set(data_list2,'Value', min( length(get(data_list2,'String')), get(data_list2,'Value') ) )
data.field = data.fields{get(data_list2,'Value')};

data_subfields_refresh;

end

%-------------------------------
function data_subfields_refresh
warning('off','sem2d_snapshot_read:FileNotFound');
data.loaded = sem2d_snapshot_read(data.field,data.time,data.dir);
warning('on','sem2d_snapshot_read:FileNotFound');
if isempty(data.loaded)
  my_error('Data not found');
  return
end
if isstruct(data.loaded)
  data.subfields = fieldnames(data.loaded);
else
  data.subfields = {''};
end
set(data_list3,'String',data.subfields);
set(data_list3,'Value', min( length(get(data_list3,'String')), get(data_list3,'Value') ) )
data.subfield = data.subfields{get(data_list3,'Value')};
end

%-------------------------------
function color_menu_callback(source,eventdata)
col = get(color_menu,'Value');
switch col
  case 1, 
    set(color_min,'Enable','off');
    set(color_max,'Enable','off');
    data.scale = 'minmax';
  case 2, 
    set(color_min,'Enable','off');
    set(color_max,'Enable','off');
    data.scale = 'maxabs';
  case 3, 
    set(color_min,'Enable','on');
    set(color_max,'Enable','on');
    cmin = str2num( get(color_min,'String') );
    cmax = str2num( get(color_max,'String') );
    data.scale = [cmin cmax];
end
end

%-------------------------------
function color_min_callback(source,eventdata)
cmin = str2num( get(color_min,'String') );
if isempty(cmin), 
  my_error('Min must be a number'); 
  set(color_min,'String',num2str(data.scale(1)));
else
  data.scale(1) = cmin;
end
end

%-------------------------------
function color_max_callback(source,eventdata)
cmax = str2num( get(color_max,'String') );
if isempty(cmax), 
  my_error('Max must be a number'); 
  set(color_max,'String',num2str(data.scale(2)));
else
  data.scale(2) = cmax;
end
end

%-------------------------------
function plot_callback(source,eventdata)
data.fig = focus_last_figure();
set(data.fig,'Name','SEM2DPACK snapshot');
% wait message
cla
lims = axis;
hwait=text((lims(1)+lims(2))/2,(lims(3)+lims(4))/2,'... Wait while plotting ...');
set(hwait,'HorizontalAlignment','center')
drawnow
% plot
if isstruct(data.loaded)
  sem2d_snapshot_plot( getfield(data.loaded,data.subfield), data.grid,data.scale);
  fieldtxt = [data.field ' ' data.subfield];
else
  sem2d_snapshot_plot(data.loaded, data.grid,data.scale);
  fieldtxt = data.field;
end
% title
DIRLEN = 30;
dirtxt = strrep(data.dir,'_','\_');
if length(dirtxt)>DIRLEN-3, dirtxt = ['...' dirtxt(end-(DIRLEN-3)+1:end)]; end
title({['dir: ' dirtxt] ,...
       ['field: ' fieldtxt], ...
       ['snapshot : ' num2str(data.time)]})
% delete wait message
delete(hwait);
drawnow
end

%-------------------------------
% get handle of last selected figure, excluding GUI figure
% if no figure is found create a new one
function h = focus_last_figure()
hs=get(0,'Children');
hs= hs( strmatch('figure',get(hs,'Type')) ); % filter out non-figures
if length(hs)<2
  h = figure;
else
  h = hs(2);
end
figure(h);
gca;
end

%-------------------------------
function close_callback(source,eventdata)
close(f)
end

%-------------------------------
function my_error(msg)
uiwait( errordlg(msg,'SEM2D_GUI error','modal') )
end

%-------------------------------
% fix directory name (full path)
function p = fix_path(inputdir)
here=pwd; 
cd(inputdir); 
p=pwd; 
cd(here);
end

%-------------------------------
end
