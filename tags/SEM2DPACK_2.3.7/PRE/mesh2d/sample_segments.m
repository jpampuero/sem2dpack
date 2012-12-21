% SAMPLE_SEGMENTS generates points that regularly sample multiple segments of a line
%
% SYNTAX	out = sample_segments(PTS,NEL)
%
% INPUTS	PTS	points defining N segments of a line or curve, size=[N+1,2]
% 		NEL 	number of elements per segment, length=N
%
% OUTPUTS	out	sampled points, size=[sum(NEL)+1,2]
%
function out = sample_segments(PTS,NEL)

out = zeros(sum(NEL)+1,2);
out(1,:) = PTS(1,:);
n = 1;
for k=1:length(NEL),
  np = NEL(k);
  shape2 = [1/np:1/np:1]';
  shape1 = 1-shape2;
  out(n+[1:np],1) = PTS(k,1)*shape1 + PTS(k+1,1)*shape2;
  out(n+[1:np],2) = PTS(k,2)*shape1 + PTS(k+1,2)*shape2;
  n = n+np;
end
