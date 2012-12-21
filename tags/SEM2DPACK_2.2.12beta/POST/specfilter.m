%specfilter.m	Zero-phase Butterworth filter
%		via spectral domain
%
% datout = specfilter(datin,samp,flo,fhi,ord)
%	datin	Real input data, length must be 2^p
%	samp	Sampling frequency
%	flo	Low frequency cut-off (0 if no low filtering)
%	fhi	High frequency cut-off (0 if no high filtering)
%	ord	Filter order (small integer)

function datout = specfilter(datin,SAMP,FLO,FHI,ORD)

if any(size(datin)==1), datin = datin(:); end
[Nin,nsis] = size(datin);
N = 2^nextpow2(Nin);
fdat = fft(datin,N);

% build the filter in spectral domain
ord2 = 2*ORD ;
% data is real, work with the positive frequencies
f = (0:N/2).' * SAMP/N ;
% filter out low frequencies
if FLO > 0
  filtLo = 1 - 1 ./ sqrt( 1 + (f/FLO).^ord2 ) ; 
else
  filtLo = ones(N/2+1,1);
end
% filter out high frequencies
if FHI > 0
  filtHi = 1 ./ sqrt( 1 + (f/FHI).^ord2 ) ; 
else
  filtHi = ones(N/2+1,1);
end
filt = filtHi .* filtLo ;

% apply the filter
fdat(1:N/2+1,:) = fdat(1:N/2+1,:) .* filt(:,ones(nsis,1)) ;
% data is real, fill the negative frequencies
fdat(N/2+2:N,:) = conj( fdat(N/2:-1:2,:) );

datout = real(ifft(fdat));
datout = datout(1:Nin,:);
