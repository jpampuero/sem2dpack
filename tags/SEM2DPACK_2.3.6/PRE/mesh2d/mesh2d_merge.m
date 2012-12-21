% MESH2D_MERGE 	merges several meshes into a single mesh
% 		Assumes each domain has 4 boundaries,
%		all domains have same type of elements (Q4 or Q9)
%
% SYNTAX	mesh=mesh2d_merge(domains,mtable)
%
% INPUTS	domains(:)	mesh structures (see MES2D_TFI) for each of the N domains
%		mtable(4,N)	merging tables:
%		 		if mtable(i,j)>0 then 
%				  merge boundary #i of domain #j to domain #mtable(i,j)
%		 		else 
%		   	  	  merge to boundary #mtable(i,j) of the merged mesh
%
% OUTPUTS	mesh	mesh structure of the merged mesh
%
function mesh=mesh2d_merge(domain,MergeTo)

[Q,NEL] = size( [domain.enod] );
ndom = length(domain);

% init local-to-global numbering 
NNOD = 0;
for idom=1:ndom,
  nnod = size(domain(idom).coor,2);
  domain(idom).iglob = zeros(nnod,1);
end

% init external boundaries
jb=find(MergeTo<=0);
ib=unique(-MergeTo(jb));
NBNDS = length(ib);
if max(ib)~=NBNDS, error('Non consecutive boundary tags'), end
for ib=1:NBNDS,
  mesh.bnds{ib} = [];
end

% generate local-to-global numbering and external boundaries
nelast=0;
for idom=1:ndom,
  inew=find(~domain(idom).iglob);
  domain(idom).iglob(inew) = NNOD+(1:length(inew));
  NNOD = NNOD+length(inew);
  nb = length( domain(idom).bnds );
  nel = size( domain(idom).enod, 2);
  for ib=1:nb,
    jdom=MergeTo(ib,idom);
    if jdom > idom
      jb=find(MergeTo(:,jdom)==idom);
      inod=GetBoundaryNodes(domain(idom),ib,Q);
      jnod=GetBoundaryNodes(domain(jdom),jb,Q);
      domain(jdom).iglob(jnod) = domain(idom).iglob(inod);
    elseif jdom <= 0
      jb = -jdom;
      newdata = domain(idom).bnds{ib};
      newdata(1,:) = nelast+newdata(1,:);
      mesh.bnds{jb} = [ mesh.bnds{jb} newdata ];
    end 
  end
  nelast = nelast+nel;
end

% generate global coordinates
mesh.coor = zeros(2,NNOD);
mesh.enod = zeros(Q,NEL);
mesh.etag = zeros(1,NEL);
nelast=0;
for idom=1:ndom,
  mesh.coor(:,domain(idom).iglob) = domain(idom).coor;
  nel = size( domain(idom).enod, 2);
  ie=nelast+(1:nel);
  mesh.enod(:,ie) = domain(idom).iglob( domain(idom).enod );
  mesh.etag(ie) = domain(idom).etag;
  nelast = nelast+nel;
end

%------------
% output is lexically sorted in polar coordinates (r,theta)
function inod = GetBoundaryNodes(domain,ib,q)

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
coor = domain.coor(1,inod) +1i*domain.coor(2,inod);
[coor,isor] = sort(coor);
inod=inod(isor);
