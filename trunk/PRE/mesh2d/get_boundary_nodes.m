function [coor,inod] = get_boundary_nodes(domain,ib,q)

if ~exist('q','var'), q=4; end
switch q
case 9,
  FaceNodes = [ [1 5 2]; [2 6 3]; [3 7 4]; [4 8 1] ]';
  qb = 3;
case 4
  FaceNodes = [ [1 2]; [2 3]; [3 4]; [4 1] ]';
  qb = 2;
otherwise
  error('Only Q4 and Q9 are implemented')
end

bdata=domain.bnds{ib};
nel = size(bdata,2);
inod = zeros(qb,nel);
for ie=1:nel,
  inod(:,ie) = domain.enod( FaceNodes(:,bdata(2,ie)), bdata(1,ie) );
end
inod = unique(inod(1:numel(inod)));
% lexicographic sort in cartesian coordinates :
[coor,isor] = sortrows(domain.coor(:,inod)');
coor = coor';
inod=inod(isor);
