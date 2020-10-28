function [it, dt, t, eqnum, isdynamic, isswitch, isEQ]=sem2d_read_rsf_timeinfo(file)
% read time step information for rate state fault with adaptive time stepping
% it: time step id
% dt: time step (s)
%  t: time (s)
% eqnum: number of earthquake
% isdynamic: if the step is a dynamic event
% isswitch: if there's a switch betweent static and dynamic at this step
%

fid       = fopen(file);
C         = textscan(fid,'%d%f%f%d%s%s%s');
it        = cell2mat(C(1));
dt        = cell2mat(C(2));
t         = cell2mat(C(3));
eqnum     = cell2mat(C(4));
isdynamic = cellfun(@(x) x=='T', C{5});
isswitch  = cellfun(@(x) x=='T', C{6});
isEQ      = cellfun(@(x) x=='T', C{7});
end
