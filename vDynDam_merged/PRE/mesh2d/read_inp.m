% READ_INP	Reads a 2D mesh in Abaqus' INP format, as output by the CUBIT mesh generation software
%		(http://cubit.sandia.gov)
%               Assumes 4-node quad elements
% 
% SYNTAX	m = read_inp(file)
%
% INPUTS	file	file name
%
% OUTPUTS	m	mesh structure in MESH2D format (see MESH2D_TFI)
%
function m = read_inp(file)

stream=fopen(file,'r');

line = search_pattern('*NODE, NSET=ALLNODES',stream);
data = fscanf(stream,'%u , %f , %f',[3 inf]);
m.coor = zeros(2,size(data,2));
m.coor(:, data(1,:) ) = data(2:3,:);

% we assume here 4-node quad elements ("S4R")
Q = 4;
line = search_pattern('*ELEMENT, TYPE=S4R',stream);
data = fscanf(stream,'%u , %u , %u , %u , %u',[1+Q inf]);
m.enod = zeros(Q,size(data,2));
m.enod(:, data(1,:) ) = data(2:Q+1,:);
% to do: input material tags
m.etag = ones(1,size(data,2));

nelset=0;
nsurf = 0;
[line,kw] = search_next_keyword(stream);
quit = 0;
m.bnds = [];
while ~quit
  switch kw
    case 'ELSET'
      nelset = nelset+1;
      elset(nelset) = read_elset(stream,line);  
      line = fgetl(stream);
    case 'SURFACE'
      nsurf=nsurf+1;
      [m.bnds{nsurf},m.bnames{nsurf},line] = read_surface(stream,line,elset);
    otherwise
      quit = 1;
  end
  [line,kw] = search_next_keyword(stream,line);
end



%---------
function line = search_pattern(pattern,stream)

line = fgetl(stream);
while isempty(strfind(line,pattern))
  line = fgetl(stream);
end

%---------
function [line,kw] = search_next_keyword(stream,line)

if ~exist('line','var'), line =[]; end
if isempty(line), line = fgetl(stream); end
kw = regexp(line,'^\*[A-Z]+','match');
while isempty( kw )
  line = fgetl(stream);
  kw = regexp(line,'^\*[A-Z]+','match');
end
kw = kw{1}(2:end);

%---------
function [elset,line] = read_elset(stream,line)

elset.name = sscanf(line,'*ELSET, ELSET=%s');
if isempty(elset.name), line = search_pattern('*ELSET',stream); end
elset.elem = fscanf(stream,'%u , ');

%---------
function [btab,bname,line] = read_surface(stream,line,elset)

bname = sscanf(line,'*SURFACE, NAME=%s');
btab = [];
if isempty(bname), line = search_pattern('*SURFACE',stream); end
line = fgetl(stream);
while line(1) ~= '*'
  stuff = textscan(line,'%s','Delimiter',',');
  elset_name = char(stuff{1}(1));
  edge_index = sscanf(char(stuff{1}(2)), 'E%u'); 
  n = strmatch(elset_name,{elset.name});
  nelem = length(elset(n).elem);
  btab = [btab [elset(n).elem(:)' ; ...
                repmat(edge_index,1,nelem)] ]; 
  line = fgetl(stream);
end
