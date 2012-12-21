% SEM2D_PLOT_GRID plots a spectral element grid
%
% SYNTAX	h=sem2d_plot_grid(g)
%
% INPUT		g	spectral element grid structure (see SEM2D_READ_SPECGRID)
%
% OUTPUT	h	handle to patch object
%
% Jean-Paul Ampuero - ampuero@gps.caltech.edu - Dec 12 2008

function h=sem2d_plot_grid(g);

p = g.ngll-1;
faces = zeros(g.nelem,4*p);
for e=1:g.nelem,
  faces(e,:) = [g.ibool(1:p,1,e)' g.ibool(p+1,1:p,e) g.ibool(p+1:-1:2,p+1,e)' g.ibool(1,p+1:-1:2,e)];
end
h=patch('Vertices',g.coord,'Faces',faces,'FaceColor','none');
axis equal
xlabel('X')
ylabel('Z')
