% READFTQ reads a Q4 mesh (4-node quadrangles) from an FTQ file
%
% SYNTAX	ftq = readftq(name)
%
% INPUT		name	prefix of the FTQ file name
%
% OUTPUT	ftq  	structure containing:
%        	ftq.npgeo		number of nodes
%        	ftq.nelem		number of elements
%        	ftq.etags(nelem)	domain tag for each element
%        	ftq.knods(nelem,4)	node indices for each element (counterclockwise)
%        	ftq.ktags(npgeo)	boundary tags for each node
%        	ftq.coorg(npgeo,2)	node coordinates
%
%  See also PLOTFTQ
%
%  Jean-Paul Ampuero jampuero@princeton.edu  May 15 2003

function ftq = readftq(name)

name = deblank(name);
if isempty(regexp(name,'\.ftq$')), name = strcat(name,'.ftq'); end
fid = fopen(name);
header = fscanf(fid,'%i',[4 1]);
ntri = header(3);
if ntri >0, error( sprintf('There are %n triangles in your mesh!',ntri) ), end

ftq.npgeo = header(1);
ftq.nelem = header(2);

elems = fscanf(fid,'%i',[6 ftq.nelem])'; % = "4" node1 node2 node3 node4 tag
nodes = fscanf(fid,'%f',[3 ftq.npgeo])'; % = X Y tag
fclose(fid);

ftq.knods = elems(:,2:5);
ftq.etags = elems(:,6);
ftq.coorg = nodes(:,1:2);
ftq.ktags = nodes(:,3)';
