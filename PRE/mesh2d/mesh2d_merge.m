% MESH2D_MERGE 	merges several meshes into a single mesh
% 		Assumes all domains have same type of elements (Q4 or Q9)
%
% SYNTAX	mesh = mesh2d_merge(domains,bnd_mtab,nod_mtab)
%
% INPUTS	domains(:)	mesh structures (see MES2D_TFI) for each of the N domains
%		bnd_mtab(:,:)	boundary merging table:
% 		 		  if bnd_mtab(i,j)>0 : merge boundary #i of domain #j to domain #bnd_mtab(i,j)
% 		 		  if bnd_mtab(i,j)<0 : merge to boundary -#bnd_mtab(i,j) of the merged mesh
%		 		  if bnd_mtab(i,j)=0 : do nothing, boundary #i of domain #j is empty
%				Merged boundaries must have same number of nodes.
%               nod_mtab(:,4)   node merging table (optional):
%                               merge node #nod_mtab(i,2) of domain #nod_mtab(i,1)
%                               with  node #nod_mtab(i,4) of domain #nod_mtab(i,3)
%
% OUTPUTS	mesh	mesh structure of the merged mesh
%
% NOTE		If the coordinates of matching node are slightly different 
%		the coordinates from the domain with lowest index prevail.
%
function mesh=mesh2d_merge(domain,MergeTo,NodesMergeTo)

if ~exist('NodesMergeTo','var'), NodesMergeTo=[]; end
% make sure the domain with lowest index is listed first
for k=1:size(NodesMergeTo,1), 
  if NodesMergeTo(k,1)>NodesMergeTo(k,3)
    NodesMergeTo(k,:) = NodesMergeTo(k,[3 4 1 2]) ;
  end
end

[Q,NEL] = size( [domain.enod] );
ndom = length(domain);

% init local-to-global numbering 
NNOD = 0;
for idom=1:ndom,
  nnod = size(domain(idom).coor,2);
  domain(idom).iglob = zeros(nnod,1);
end

% init external boundaries
jb=find(MergeTo<0);
ib=unique(-MergeTo(jb));
NBNDS = max(ib);
for ib=1:NBNDS,
  bnds{ib} = [];  
  % we don't initialize mesh.bnds at this point
  % because bnds must be the last item in the structure
end

% generate local-to-global numbering and external boundaries
nelast=0;
for idom=1:ndom,
  
 % nodes that have not yet been processed are given a new number
  inew=find(~domain(idom).iglob);
  domain(idom).iglob(inew) = NNOD+(1:length(inew));
  NNOD = NNOD+length(inew);

 % merge nodes
  for k=1:size(NodesMergeTo,1),
    if idom ~= NodesMergeTo(k,1), continue; end
    inod = NodesMergeTo(k,2);
    jdom = NodesMergeTo(k,3);
    jnod = NodesMergeTo(k,4);
    domain(jdom).iglob(jnod) = domain(idom).iglob(inod);
  end

 % merge the nodes of each boundary of the current domain
 % with matching nodes in the matching domain
 % or add boundary elements to a boundary of the global domain
  nb = length( domain(idom).bnds );
  nel = size( domain(idom).enod, 2);
  for ib=1:nb,
    if isempty( domain(idom).bnds{ib} ), continue; end
    jdom=MergeTo(ib,idom);
   % merge nodes to matching boundary of a different domain
    if jdom > idom 
      jb=find(MergeTo(:,jdom)==idom);
      inod=GetBoundaryNodes(domain(idom),ib,Q);
      jnod=GetBoundaryNodes(domain(jdom),jb,Q);
      if (length(inod)~=length(jnod)), 
        error(['Boundary ' num2str(ib) ' of domain ' num2srt(idom) ...
          ' and boundary ' num2str(jb) ' of domain ' num2srt(jdom) ...
          ' have different number of nodes']), end 
      domain(jdom).iglob(jnod) = domain(idom).iglob(inod);
   % add boundary elements to a boundary of the global domain
    elseif jdom < 0 
      jb = -jdom;
      newdata = domain(idom).bnds{ib};
      newdata(1,:) = nelast+newdata(1,:);
      bnds{jb} = [ bnds{jb} newdata ];
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
mesh.bnds = bnds;

%------------
% output is lexically sorted 
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
% lexicographic sort in cartesian coordinates :
%[coor,isor] = sortrows(domain.coor(:,inod)');
dmax1 = max(domain.coor(1,inod)) - min(domain.coor(1,inod));
dmax2 = max(domain.coor(2,inod)) - min(domain.coor(2,inod));

if dmax1>dmax2
    [~, isor] = sort(domain.coor(1,inod));
else
    [~, isor] = sort(domain.coor(2,inod));
end

%coor = coor';
% lexicographic sort in polar coordinates (r,theta) :
%coor = domain.coor(1,inod) +1i*domain.coor(2,inod);
%[coor,isor] = sort(coor);
inod=inod(isor);
