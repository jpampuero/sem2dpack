% [Trup,Tpz] = PlotRuptureFront(V,Veps,D,Dc,X,DT)
function [Trup,Tpz] = PlotRuptureFront(V,Veps,D,Dc,X,DT)

NX = length(X);
[N1,N2] = size(V);
if N1==NX
  NT=N2;
elseif N2==NX
  V = V';
  D = D';
  NT=N1;
else
  error('Size mismatch (V/D/X)')
end

Trup = repmat(NT,NX,1);
Tpz  = repmat(NT,NX,1);

% force V(:,1)=0 at t=0
%if nnz(V(:,1))
%  V = [zeros(NX,1) V];
%  D = [zeros(NX,1) D];
%end

for k=1:NX,

 % rupture front
  m = find( V(k,:)>Veps );
  if ~isempty(m)
  if m>1
    m = m(1)-1;
    Trup(k) = m-1 + (Veps-V(k,m))/(V(k,m+1)-V(k,m));
  else
    Trup(k) = 0;
  end
  end

 % end of process zone
  m = find( D(k,:)<Dc);
  m = m(end);
  if m<NT
    Tpz(k) = m-1 + (Dc-D(k,m))/(D(k,m+1)-D(k,m));
  end

end

Trup = Trup*DT;
Tpz = Tpz*DT;

if nargout==0
  %plot(X,Tpz,X,Trup)
  area( X,Tpz,'FaceColor','b')
  hold on; area( X,Trup,'FaceColor','w'); hold off
  xlabel('X')
  ylabel('T')
end
