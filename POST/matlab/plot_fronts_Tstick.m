% PLOT_FRONTS_TSTICK space-time plot of rupture front and process zone tail 
% using T "stick"
%
% SYNTAX	[Trup,Tpz] = plot_fronts(front_mu, front_mustick, mu_s, D, Dc, X, DT)
%
% INPUTS	
%		front_mu(:,:)	space-time friction coefficient data
%		front_mustick(:,:)	space-time "stick" friction coefficient
%		(tau_stick/siigma_n)
%		mu_s	static friction coefficient
%		D(:,:)	space-time slip data
%		Dc	slip threshold to define the tail of process zone
%		X(:)	positions of slip and velocity data
%		DT	timestep of slip and velocity data
%		
% OUTPUTS	
%       Trup(:) rupture times (for which tau=tau_s interpolated from tau stick)
%		Tpz(:)  time of tail of process zone (for which D=Dc) 
%			at each location, with quadratic time interpolation
%		Plot  	a shaded area in space-time indicating the head and 
%			tail of the rupture front
%
function [Trup,Tpz] = plot_fronts_Tstick(mu, mu_st, mu_s, D, Dc, X, DT)

NX = length(X);
[N1,N2] = size(mu);
if N1==NX
  NT=N2;
elseif N2==NX
  mu = mu';
  D = D';
  NT=N1;
else
  error('Size mismatch (V/D/X)')
end

Trup = repmat(NT,NX,1);
Tpz  = repmat(NT,NX,1);

for k=1:NX
    
      % rupture front (interpolate from mu stick)
      m = find(mu_st(k,:) > mu_s, 1, 'first');
      if isempty(m)
          m = find(mu(k,:) == max(mu(k,:))); 
          if size(m,2)~=1 
                Trup(k) = 0;
          else
                Trup(k) = m;
          end
      else
          Trup(k) = interp1(mu_st(k,m-1:m), m-1:m, mu_s, 'linear', 'extrap'); % assumes m=1 is t=0
          
      end
    
    
 % end of process zone
  m = find( D(k,:) < Dc);
  m = m(end);
  if m<NT
    Tpz(k) = m-1 + (Dc-D(k,m))/(D(k,m+1)-D(k,m)); % assumes m=1 is t=0
  end

end

Trup = Trup*DT;
Tpz = Tpz*DT;

if nargout==0
  plot(X,Tpz,X,Trup)
  legend('Rupture front', 'Rupture tail');
  xlabel('X')
  ylabel('T')
end
