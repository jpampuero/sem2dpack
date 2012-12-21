% PLOTFTQ plots a Q4 mesh from FTQ data
%
% SYNTAX	plotftq(ftqdata)  
%
% INPUT		ftqdata	an FTQ structure, as defined in READFTQ
%        	mode	plot mode, can be either of:
%			  1 = plots only the mesh (default)
%      			  2 = colors each element according to its domain tag
%
%  Jean-Paul Ampuero jampuero@princeton.edu  May 15 2003

function h=plotftq(ftq,mode)

if ~exist('mode','var'), mode=1; end

clf
if mode==1
  h=patch('Vertices',ftq.coorg,'Faces',ftq.knods,'FaceColor','none');
else
  colormap('cool');
%  set(gca,'DefaultPatchEdgeColor','none');
  h=patch('Vertices',ftq.coorg,'Faces',ftq.knods,'FaceVertexCData',ftq.etags,'FaceColor','flat');
end
axis equal
xlabel('X'); ylabel('Y');
