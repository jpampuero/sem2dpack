% SEM2D_EXTRACT_POINT extracts field value at an arbitrary point, by SEM interpolation
%
% SYNTAX	F = sem2d_extract_points(field,g,xz)
%
% INPUT		field(:,:,:)	field in local storage
%		g	grid structure
%		xz(2)	coordinates of the target point	
%
% OUTPUT	F	value of field at point xz
%			interpolated by SEM basis functions
%
function F = sem2d_extract_point(field,g,xz)

% find global index of nearest node
dist2 = (g.coord(:,1)-xz(1)).^2 + (g.coord(:,2)-xz(2)).^2 ;
[tmp,k] = min(dist2);

% get local indices of nearest node(s)
% the nearest node might belong to several elements
ijes = find(g.ibool(:)==k);
[is,js,es] = ind2sub( size(g.ibool), ijes);

% find local coordinates (xi,eta,e)
for n=1:length(es),
  [xi,eta,xz_new,istat] = fe_find_point( xz, g, es(n), g.x(is(n)), g.x(js(n)) );
  if istat==0, break, end	% exit if found node in this element 
end
e = es(n);

% compute interpolation tables
interp = se_init_interpol(xi,eta,g.x,g.ngll);

% interpolate field
F = field(:,:,e);
F = interp(:)' * F(:);



%----------
function [xi,eta,coord_new,istatus] = fe_find_point( coord, g, e, xi0,eta0 )

NTRIAL = 100;
TINY_XABS = 1e-9;

coord = coord(:);
coorg = g.coorg(g.knods(e,:),:)';
x = [xi0; eta0];
tolx = TINY_XABS;
tolf = TINY_XABS*sqrt( fe_element_area(coorg) );

for n=1:NTRIAL,

  shap = q4_elem('shape',x);
  fvec = coorg*shap - coord;

  dshap = q4_elem('dershape',x);
  fjac = coorg * dshap;

  if all(abs(fvec) <= tolf), break; end
  p = - fjac \ fvec;
  x = x + p;
  if all(abs(p) <= tolx), break; end

end

%if n>=NTRIAL, error('fE_find_point: did not converge'); end
if n>=NTRIAL, istatus=2; end
if all( abs(x)< 1.d0+TINY_XABS )
  istatus = 0;
else
  istatus = 1;
end
xi=x(1);
eta=x(2);
shap = q4_elem('shape',x);
coord_new = coorg * shap;


%----------------------
function area = fe_element_area(coorg)

a = (coorg(1,2)-coorg(1,1))^2+(coorg(2,2)-coorg(2,1))^2;
b = (coorg(1,3)-coorg(1,2))^2+(coorg(2,3)-coorg(2,2))^2;
c = (coorg(1,4)-coorg(1,3))^2+(coorg(2,4)-coorg(2,3))^2;
d = (coorg(1,1)-coorg(1,4))^2+(coorg(2,1)-coorg(2,4))^2;
p = (coorg(1,1)-coorg(1,3))^2+(coorg(2,1)-coorg(2,3))^2;
q = (coorg(1,2)-coorg(1,4))^2+(coorg(2,2)-coorg(2,4))^2;
area = 0.25*sqrt( 4*p*q - (b+d-a-c)^2 );


%----------------------
function out = q4_elem(mode,x)

s=x(1);
t=x(2);

switch lower(mode)
  case 'shape'

  sp = s + 1;
  sm = s - 1;
  tp = t + 1;
  tm = t - 1;

  out(1) = sm * tm;
  out(2) = - sp * tm;
  out(3) = sp * tp;
  out(4) = - sm * tp;
  
  out = out(:)/4;

  case 'dershape'

  sp = s + 1;
  sm = s - 1;
  tp = t + 1;
  tm = t - 1;

  out(1,1) = tm;
  out(2,1) = - tm;
  out(3,1) =  tp;
  out(4,1) = - tp;

  out(1,2) = sm;
  out(2,2) = - sp;
  out(3,2) =  sp;
  out(4,2) = - sm;

  out = out/4;

end


%----------------------
function interp = se_init_interpol(xi,eta,xgll,ngll)

k=0;
for j=1:ngll,
  fj = hgll(eta,xgll(j),ngll);
  for i=1:ngll,
    k=k+1;
    fi = hgll(xi,xgll(i),ngll);
    interp(k) = fi*fj;
  end
end

%=====================================================================
% The following functions are translated from gll.f90

  function HGLL = hgll(Z,ZGLL,NZ)

  EPS = 1e-5;
  DZ = Z - ZGLL;
  if abs(DZ) < EPS 
   HGLL = 1;
   return
  end
  N = NZ - 1;
  ALFAN = N*(N+1);
  HGLL = - (1-Z*Z) * pndleg(Z,N) / ( ALFAN*pnleg(ZGLL,N)*(Z-ZGLL) );

%---------------------------------------------------------------------

  function PNDLEG = pndleg(Z,N)

  P1   = 1;
  P2   = Z;
  P1D  = 0;
  P2D  = 1;
  P3D  = 1;
  for K = 1:N-1,
   FK  = K;
   P3  = ((2*FK+1)*Z*P2 - FK*P1)/(FK+1) ;
   P3D = ((2*FK+1)*P2 + (2*FK+1)*Z*P2D - FK*P1D) /(FK+1);
   P1  = P2;
   P2  = P3;
   P1D = P2D;
   P2D = P3D;
  end
  PNDLEG = P3D;

%---------------------------------------------------------------------

  function PNLEG = pnleg(Z,N)

  P1   = 1;
  P2   = Z;
  P3   = P2;
  for K = 1:N-1,
   FK  = K;
   P3  = ((2*FK+1)*Z*P2 - FK*P1)/(FK+1);
   P1  = P2;
   P2  = P3;
  end
  PNLEG = P3;
