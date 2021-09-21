function [fas,f,v2] = FAS(t,v, vtaperstart, vtaperdur)       
%
% WHAT SEE_SV_SPECT73 DOES:  This routine inputs a slip-velocity time series
% in sq.  It deletes the leading zero-velocity part of the time series (the
% start of the truncated time series is set to 0.015s before the time when
% slip velocity > 0.01 m/s). The time series is differentiated to acceleration
% and a cosine taper is applied starting at time ataperstart and ending at 
% time (ataperstart+ataperdur).  The acceleration time series is padded 
% with zeros to extend to time (vtaperstart+vtaperdur).  The time series is
% integrated to velocity and cosine tapered from vtaperstart to
% (vtaperstart+vtaperdur).  Then ffts are taken of the tapered velocity
% time series, and plots are made. This procedure of double tapering
% extends the velocity time series with a long duration constant level,
% which is then tapered with a (hopefully) long duration taper. 

% t = tSR;
% v = U1SS(:,10);
% t = V1RAD_case0(:,1);
% v = V1RAD_case0(:,10);

% t = tsamp;
% v = vpsampSR;
% taperlengthpercent = 5;
% vtaperstart = tsamp(end)*(100-taperlengthpercent)/100;
% vtaperdur = tsamp(end)-vtaperstart;

Dt = t(2)-t(1);
Fs = 1/Dt;
v2 = v;


ataperstart = vtaperstart;
ataperdur = vtaperdur;
minfreq = 0.005;

dt = t(2)-t(1);

i1 = 1;%find(v>1e-4, 1 ); 
t1 = t(i1); % time of the first sample with nonzero v
tend = t(length(t)); % time of the last sample in the time series
taperend = vtaperstart + vtaperdur;
texcess = (t1+taperend) - tend; 

if texcess > 0 
    error(['see_sv_spect1: your specified taperstart and taperdur must be shortened by ' ...
        num2str(texcess) ' s'])
end


nflat = floor(vtaperstart/dt);  % number of samples which will not be tapered
ntaper = floor ( vtaperdur / dt) ; % number of samples to taper
itaper = (nflat+1):(nflat+ntaper); % samples to taper


% remove initial zero slip values from v
v(1:(i1-1)) = [];

% prepend zero slip velocity
v = [0; v];

t = dt*( 0:(length(v)-1) ); 

ataperend = ataperstart+ataperdur; 

% differentiate to acceleration
a = diff(v) / dt ; 
na = length(a); 

% discard any part of the accelerogram after ataperend because the edge effects
% enter the problem at the in-plane receivers
% nt6 = round(ataperend / dt) ; 
% if nt6 < na 
%     a( (nt6+1):na)  = [];
%     t( (nt6+1):na)  = [];
%     na = nt6; 
% else
%     error('see_sv_spect73: nt6 >= na.  This option not yet coded.')
% end

% now the acceleration has na points, not exceeding ataperend in duration
% taper the last ataperdur s of the acceleration time series
nataper = round(ataperdur/dt);
% calculate the taper for acceleration 
taper = 0.5 * ( 1 + cos( (1:nataper) * pi / nataper)); % taper weights
jtaper = (na-nataper+1):na ; % indices to taper
auntapered = a; 
a(jtaper)= a(jtaper) .* taper' ; % apply taper to acceleration
% a(na) is now zero.  Add more zeros
ntapertotal = nflat+ntaper; 
kk = (na+1):ntapertotal ; 
a( kk ) = 0; 
t( kk ) = (kk-1)*dt; 


% total length of time series is now nflat+ntaper, 
% integrate to velocity
v = cumsum(a) * dt; 
t = dt*( 0:(length(v)-1) ); 

% apply taper to velocity
taper = 0.5 * ( 1 + cos( (1:ntaper) * pi / ntaper)); % taper weights
v(itaper) = v(itaper) .* taper'; % apply taper

v = v2;


% remove the remaining points after taper
tv = v(1:(nflat+ntaper)); 
ntv = length(tv); 

% normalize
tv = tv ;%/ ( max(abs(tv)) ); 

% add zeros to make length a power of 2
nfft = 2^nextpow2(ntv); 
mfft = nfft ; % retain info for plotting

nfft = 2^nextpow2(1/minfreq/dt); % I multiply the number of time samples to increase the time 
% series length, reducing df (see below).  This will give me more samples
% in the spectrum at low frequency


tvp = zeros(nfft,1); 
tvp(1:ntv) = tv; % tvp is tv padded with trailing zeros

ft = fft(tvp, nfft); % fourier transform 
% ft = ft * dt ; % normalize
fas = abs(ft); % magnitude of the fourier transform

nf = nfft/2 + 1; % number of frequencies from 0 to nyquist
df = 1 / (nfft*dt); % frequency increment is 1/ duration of time series
f = (0:(nf-1)) * df; % frequencies corresponding to each element of ft
fas((nf+1):nfft) = []; % remove elements correponding to frequencies above the nyquist

fas = fas';

% loglog( f, fas*10^shift)
% xlim([0.07 filt]);

% return

