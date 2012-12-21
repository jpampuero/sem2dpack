% READ_DCM	Reads a 2D mesh in DCM format, the output of the EZ4U mesh generation software
%		(http://www-lacan.upc.es/ez4u.htm)
%		Assumes Q4 or Q9 elements.
% 
% SYNTAX	m = read_dcm(file)
%
% INPUTS	file	file name
%
% OUTPUTS	m	mesh structure in MESH2D format (see MESH2D_TFI)
%
function m = read_dcm(file)

stream=fopen(file,'r');
[np, nel] = readHeader(stream);
m.coor = readNodes(stream, np);
m.enod = readConnectivities(stream, nel);
m.etag = ones(size(m.enod,2),1); % to do: element tags
m.bnds = readBoundaries(stream);
status=fclose(stream);

%---------
function [ nOfNodes, nOfElements ] = readHeader( stream )

nOf = fscanf(stream,'%d', 3);
nOfNodes=nOf(1);
nOfElements=nOf(2);

%---------
function  X  = readNodes( stream, nOfNodes )

X=zeros(2,nOfNodes);

for iNode=1:nOfNodes

    nodeId=fscanf(stream,'%d',1);
    x3=fscanf(stream,'%f',3);
    %X(iNode,:)=x3(1:2);
    X(:,nodeId)=x3(1:2);

end

%---------
function [enod, Q]  = readConnectivities( stream, nOfElements )

for iElement=1:nOfElements

    elementId=fscanf(stream,'%d',1); 
    dimension=fscanf(stream,'%d',1);
    p=fscanf(stream,'%d',1);
    nOfInnerNodes=fscanf(stream,'%d',4)';
    nOfNodes=nOfElementNodes( dimension, p, nOfInnerNodes);
    
    if iElement==1
      if (nOfNodes~=4 & nOfNodes~=9), error('Only Q4 and Q9 elements are allowed'); end
      Q = nOfNodes;
      enod = zeros(Q,nOfElements);
    end

    if(nOfNodes~=Q), error(['Element ' num2str(iElement) ' does not have ' num2str(Q) ' nodes']); end
        
    enod(:,elementId) = fscanf(stream,'%d',nOfNodes);

end

%---------
function nOfElementNodes = nOfElementNodes( dimension, p, nOfInnerNodes )

if(dimension==1)
    nOfElementNodes=p+1;
elseif(dimension==2)
    nOfElementNodes=nOfInnerNodes(1)+nOfInnerNodes(1)*nOfInnerNodes(2)+nOfInnerNodes(3);
end

%---------
function b = readBoundaries(stream)

b =[];
[binfo,count] = fscanf(stream,'%d',5);
while count>0
  nb = binfo(4);
  tag = binfo(5);
  tmp = fscanf(stream,'%d',[3 nb]);
  b{tag} = tmp(1:2,:);
  [binfo,count] = fscanf(stream,'%d',5); 
end
