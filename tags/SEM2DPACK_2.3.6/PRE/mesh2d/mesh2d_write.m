% MESH2D_WRITE 	writes a 2D mesh database file (*.mesh2d)
%
% SYNTAX	mesh2d_write(mesh,fname)
% 
% INPUTS	mesh	a mesh2d structure (see MESH2D_TFI)
%		fname	mesh file name 
%			The suffix .mesh2d is automatically added 
%			to the file name, if absent.
%
% OUTPUTS	fname.mesh2d	a MESH2D file. The format is documented in the
%			User's Manual of SEM2DPACK (&MESH_MESH2D input block).
%
% 

function mesh2d_write(m,fname)

fname = deblank(fname);
if isempty(regexp(fname,'\.mesh2d$')), fname = strcat(fname,'.mesh2d'); end

fid = fopen(fname,'w');
fprintf(fid,'NEL NPEL NNOD NBC\n');
fprintf(fid,'%u %u %u %u\n',size(m.enod,2),size(m.enod,1),size(m.coor,2),length(m.bnds));
fprintf(fid,'NID X Y\n');
fprintf(fid,'%u %f %f\n',[[1:size(m.coor,2)]; m.coor]);
fprintf(fid,'EID NODES TAG\n');
fmt = [repmat('%u ',1,size(m.enod,1)+2) '\n'];
fprintf(fid,fmt,[[1:size(m.enod,2)]; m.enod; m.etag]);
for k=1:length(m.bnds),
  if isempty(m.bnds{k}), continue, end
  fprintf(fid,'BCTAG NBEL\n');
  fprintf(fid,'%u %u\n',k,size(m.bnds{k},2));
  fprintf(fid,'BID EID EDGE\n');
  fprintf(fid,'%u %u %u\n',[[1:size(m.bnds{k},2)]; m.bnds{k}]);
end
fclose(fid);
