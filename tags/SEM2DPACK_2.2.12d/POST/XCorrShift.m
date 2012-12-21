% [dsamp,maxc,xcor,dsubsamp] = XCorrShift(s1,s2,win)
%
% PURPOSE:	Cross-correlation picking
%		Slides signal s1 through signal s2 to find
%		the relative delay (max normalized xcorr)
%		FFT is used for the cross-correlation, so it is optimal
%		only if input lengths are of same order.
%
% INPUT:	s1	reference signal (short) window 
%		s2	target signal (long) window    
%		win	window for search of max coherency (in samples) [] 
%
% OUTPUT:	dsamp	relative delay (in samples), 
%			ALWAYS >=0
%			0 MEANS ALIGNED
%		maxc	optimal coherency (can be negative) [-1:1]
%		xcor 	cross-correlation time series = Sxy/sqrt(Sxx*Syy)
%		dsubsamp sub-sample precision delay 
%			estimated by quadratic interpolation
%
% NOTE:		If s1 and s2 are windows of original data,
%		picked at p1 and p2 respectively,
%		say s1 = sis1(p1+1:p1+n1) 
%		and s2 = sis2(p2+1:p2+n2),
%		then the new pick sample p2_new=p2+dsamp
%
function [dsamp,maxc,xcor,dsubsamp] = XCorrShift(s1,s2,win,only_positive)

n1=length(s1);
n2=length(s2);
if n1>n2, error('Reference signal S1 cannot be longer than target S2'), end
nx = n2-n1+1; % length of non-spurious part of xcor (no wraparound) 
if nargin<3, win=[]; end
if isempty(win)
  win=[1:nx]; 
else
  win=[max(1,win(1)):min(nx,win(end))];
end
if nargin<4, only_positive=0; end

% cross-correlation
nfft = 2^nextpow2(n2+n1);
f1=fft(s1,nfft);
f2=fft(s2,nfft);
xcor=real(ifft( conj(f1).*f2 ));
xcor=xcor(1:nx);

% scale
s2s2 = s2.*s2;
scal=zeros(nx,1);
scal(1) = sum(s2s2(1:n1));
for k=1:nx-1,
  scal(k+1) = scal(k)+ s2s2(n1+k)-s2s2(k);
end
scal = sqrt(scal)*norm(s1);
xcor=xcor./scal ;

% optimal lag
if only_positive
  [maxc,dsamp]=max(xcor(win));
else
  [maxc,dsamp]=max(abs(xcor(win)));
end
dsamp = dsamp +win(1)-1;
maxc=xcor(dsamp);

if nargout>3
if dsamp>1 & dsamp<nx-1
  dsubsamp = dsamp-0.5*( xcor(dsamp+1)-xcor(dsamp-1) ) ...
             /( xcor(dsamp-1)-2*xcor(dsamp)+xcor(dsamp+1) );
  dsubsamp = dsubsamp-1;
else
  dsubsamp = dsamp-1;
end
end

% lag --> delay
dsamp= dsamp-1;
