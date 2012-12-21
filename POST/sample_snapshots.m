% SAMPLE_SNAPSHOTS example of visualization of snapshot outputs from SEM2DPACK
% You must adjust the parameters FNAME and ISNAP in this script

% select a field : 
FNAME = 'vx';
%			PSV			SH
% 	displacement	dx,dz			dy
%	velocity	vx,vz			vy
%	acceleration	ax,az			ay
%	strain		e11,e22,e12 		e13,e23
%	stress		s11,s22,s33,s12		s13,s23

% select snapshot indices:
ISNAP = [1:20];

grid = sem2d_read_specgrid();

for is=ISNAP,
  field = sem2d_snapshot_read(FNAME,is);
  if isempty(field), break, end
  clf
  sem2d_snapshot_plot(field,grid);
  title( [FNAME ' snapshot # ' num2str(is)])
  pause
end
