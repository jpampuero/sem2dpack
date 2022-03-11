% SEM2D_EXTRACT_LINE extract values of a field along a vertical or horizontal line
%	 	from a box mesh (rectangular cartesian grid), and sorts them by position
%
% SYNTAX	[XZ,F,xz_out] = sem2d_extract_line(field,grid,xz,hv)
%
% INPUT		field 	data field in local storage format, see sem2d_snapshot_read
%		grid 	spectral element grid structure, see sem2d_read_specgrid 
%		xz	target value of the line position (x-position if vertical, 
%			z-position if horizontal). The output is the nearest line of nodes.
%		hv	line orientation: 'h' for horizontal, 'v' for vertical line
%
% OUTPUT	XZ	position of nodes in the line (x or z, as in input)
%		F	fields values at nodes
%		xz_out	value of xz possibly relocated to nearest node
%
% EXAMPLES	Plot damage variable along two horizontal lines, right above and right below
%		a horizontal fault located at x=0:
%			g=sem2d_read_specgrid;
%  			d=sem2d_snapshot_read('dmg',10);
%			[Xu,Fu] = sem2d_extract_line(d.alpha,g,0.0001,'h');
%			[Xd,Fd] = sem2d_extract_line(d.alpha,g,-0.0001,'h');
%			plot(Xu,Fu, Xd,Fd)
%
%		Plot damage along a vertical line at x=15, as a function of distance to the fault
%		separating the sections below and above the fault:
%			[Z,F] = sem2d_extract_line(d.alpha,g,15,'v');
%			k = find(Z==0);
%			semilogx(-Z(1:k(1)),F(1:k(1)), Z(k(2):end),F(k(2):end))
%			legend('below','above')

function [X,F,xz_nearest] = sem2d_extract_line(field,g,xz_target,hv)

switch hv
  case 'h',
    tdim = 1;
    ndim = 2;
    xname = 'Position along fault';
  case 'v'
    tdim = 2;
    ndim = 1;
    xname = 'Position normal to fault';
  otherwise
    error('hv must be "h" or "v"')
end

% Step 1: find the local indices of nodes along the line

% mesh information, assumes square elements
xel = g.coord( g.ibool(1:g.ngll,1,1), 1);
xel = xel-xel(1);
h = xel(end);
XTOL = 1e-5*h; 	% distance tolerance

% find the nodal x (or z)nearest to xz_target
xmin = min(g.coord(:,ndim));
x1_nearest = xmin + h*floor((xz_target-xmin)/h);
[tmp,k] = min( abs(x1_nearest+xel-xz_target) );
xz_nearest = x1_nearest + xel(k);
if abs(xz_target-xz_nearest)>XTOL
%  warning on verbose
  warning('sem2d_extract_line:relocate', ...
    sprintf('Relocated the target position %f to nearest node at %f',xz_target,xz_nearest))
end

% x of element "center" node nearest to xz_target
kc = floor(g.ngll/2);
xc_nearest = x1_nearest + xel(kc);
% WARNING what happens if the line is at the edge of the mesh?

% find elements that contain the target line position
ic = g.ibool(kc,kc,:);	% element center nodes
xc = g.coord(ic,ndim);
zc = g.coord(ic,tdim);
elist = find( abs(xc-xc_nearest)<XTOL );
% WARNING if the line is in between two elements, this only takes elements on one side

% sort element list by position
[tmp,esor] = sort(zc(elist));  
elist = elist(esor);

% loop over all elements
nodes = [];
es = zeros(g.ngll,1);
for e= elist(:)',

 % find nodes that match the target x
  x = g.coord(g.ibool(:,:,e),ndim); 	% a vector
  ks = find( abs(x-xz_nearest)<XTOL );

 % assume nodes are sorted by position (lexicographic on x,z)
 % otherwise sort nodes by position:
 % z = g.coord(g.ibool(:,:,e),tdim);
 % [tmp,isor] = sort( z(ks) );
 % ks = ks(isor);
 
 % store nodes
%  [is,js] = ind2sub( [g.ngll g.ngll], ks ); % inlined below
  js = floor((ks-1)/g.ngll)+1;
  is = ks - (js-1)*g.ngll;
  es(:) = e;
  nodes = [nodes;[is,js,es]];

end

% Step 2: extract field values at the nodes
sloc = size(g.ibool);
k = sub2ind(sloc,nodes(:,1),nodes(:,2),nodes(:,3));
iglob = g.ibool(k);
X = g.coord(iglob,tdim);
F = field(k);

