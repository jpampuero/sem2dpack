addpath ~/WAVES/SEM2DPACK_2.x/POST/  % check this path in your SEM2DPACK installation directory
sample_seis

%export at this locations:
ista = [4 18 23 46];
xsta(ista)

for n=1:length(ista),
  k=ista(n);
  subplot(2,2,n)
%  plot(t*1e6, ux(:,k), t, uz(:,k))
  plot(t*1e6, uz(:,k) )
  axis([15 100 -inf inf])
  xlabel('Time (\musec)')
  cosa=[t' ux(:,k) uz(:,k)];
  save(['p' num2str(n) '.tab'],'cosa','-ascii')
end
