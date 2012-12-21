% MESH2D_PLOT plots a 2D mesh
%
% SYNTAX	mesh2d_plot(mesh,mode)
%		
% INPUT		mesh	a mesh structure (see MESH2D_TFI).
%			If mesh(:) is a vector, the meshes of several domains
%			are plotted with common color conventions for boundaries
%               mode    plot mode, can be either of:
%                         1 = plots only the mesh (default)
%                         2 = colors each element according to its domain tag
%
% OUTPUT	A mesh plot. Each domain and each boundary has a different color.
%		The colors cycle through 'brgcym': domain or boundary #3 is green, 
%		domain or boundary #7 is blue.
%
function mesh2d_plot(domains,mode)

if ~exist('mode','var'), mode=1; end
colors='brgcym';
ncol = length(colors);
for idom=1:length(domains),
  plot_domain(domains(idom),colors(mod(idom-1,ncol)+1),colors,mode);
end

%---
function plot_domain(d,ecolor,bcolors,mode)

if size(d.enod,1)==4
% Q==4
  enod = d.enod';
  np=2; % number of nodes per element edge
  bk = [[1;2] [2;3] [3;4] [4;1]]; % local indices of nodes on each element edge
else 
% Q==9
  enod = d.enod([1 5 2 6 3 7 4 8],:)';
  np=3; 
  bk = [[1;5;2] [2;6;3] [3;7;4] [4;8;1]];
end

% plot elements
if mode==1
  patch('Vertices',d.coor','Faces',enod, 'FaceColor','none', 'EdgeColor',ecolor);
else
  colormap('cool');
  patch('Vertices',d.coor','Faces',enod, ...
        'FaceVertexCData',d.etag(:),'FaceColor','flat', ...
        'EdgeColor',ecolor);
end

% plot boundaries
ncol = length(bcolors);
for ib=1:length(d.bnds),
  if isempty(d.bnds{ib}), continue, end
  neb = length(d.bnds{ib});
  x = zeros(np,neb);
  y = zeros(np,neb);
  for eb=1:neb,
    e=d.bnds{ib}(1,eb);	% element id
    s=d.bnds{ib}(2,eb);	% edge id
    k = d.enod(bk(:,s),e); % global ids of edge nodes
    x(:,eb) = d.coor(1,k);
    y(:,eb) = d.coor(2,k);
  end
  line(x,y,'Color',bcolors(mod(ib-1,ncol)+1),'LineWidth',2)
end  
