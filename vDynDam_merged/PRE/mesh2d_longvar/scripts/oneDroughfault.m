% fractal fault profile
% November 13, 2012

% fault length L, grid spacing h=L/N
% periodicity is assumed for profile, so y(L/2)=y(-L/2)

% NODE: fault will have ODD number N+1 grid points
% x = [-L/2:h:L/2] (length L, grid spacing h=L/N)
% Fourier method will give N points with replica length L

function [x,y,n,s,res] = oneDroughfault(h,N,alpha,seed)

% h = 30; % grid spacing, nominally h=0.1 km for r=1
% N = 600;
% alpha = 10^(-3); % roughness
cutoff = true; Lmin = 5*h; % minimum wavelength (nominally 20*h)

L = N*h; % profile length
x = [-L/2:h:L/2]; % coordinates, including two endpoints

shift_profile = true; % shift profile in y direction so that yend=0

% wavenumber
  
kN = pi/h; % Nyquist wavenumber
k = zeros(1,N);
k(1:N/2+1) = 2*[0:N/2]/N;
k(N:-1:N/2+2) = -k(2:N/2);
k = k*kN;

% white noise, unit normal distribution

% seed = 7; %if you want that to be always the same
if isempty(seed); 
    s = RandStream('mt19937ar','Seed','shuffle');
else
    s = RandStream.create('mt19937ar','seed',seed);
end
% 
RandStream.setGlobalStream(s);

y = randn([1,N]);

% scale so PSD has unit amplitude

y = y*sqrt(N/L);

% FFT

Y = fft(y)*h;

% calculate PSD and check for unit amplitude

PSDy = abs(Y).^2/L;
% disp(['PSD = ' num2str(mean(PSDy))])

% multiply Y by square-root of desired PSD

PSDy_exact = (2*pi)^3*alpha^2*k.^(-3);
Y = Y.*sqrt(PSDy_exact);

% remove k=0 component

Y(1) = 0;

% add short wavelength (high wavenumber) cutoff

if cutoff
  kmax = 2*pi/Lmin;
  I = find(abs(k)>kmax); 
  Y(I) = 0;
end

% inverse FFT

y = ifft(Y)/h;
y = real(y); % just to clean up

% calculate slope

M = 1i*k.*Y;
M(N/2+1) = 0; % no Nyquist for spectral differentiation
m = ifft(M)/h;
m = real(m);
  
% convert slope into unit normals

n = zeros(N,2);
n(:,1) = -m./sqrt(1+m.^2);
n(:,2) = 1./sqrt(1+m.^2);

% check alpha

alpha_check = sqrt(h*sum(y.^2)/L)/L; res.alpha=alpha_check;
res.ratio = alpha/alpha_check;
res.per = 100*abs(alpha-alpha_check)/alpha;
% disp(['input alpha=',num2str(alpha),' calculated alpha=',num2str(alpha_check)])
% disp(['ratio=',num2str(ratio), ', difference %', num2str(per*100)])

% only return non-negative portion of spectrum

k = k(1:N/2+1); k(N/2+1) = kN;
Y = Y(1:N/2+1);
M = M(1:N/2+1);
PSDy = abs(Y).^2/L;
PSDy_exact = PSDy_exact(1:N/2+1);

%return

% repeat first point, shift if needed

y = [y y(1)];
n = [n; n(1,:)];

xend = [-L/2 L/2]; % endpoints of fault
yend = y(1); if shift_profile, y = y-yend; yend = 0; end

% plot unrotated profile (with vertical exaggeration)

% figure(1),clf
% 
% plot(x,y,'b+');
% % daspect([1 1e-5 1]);
% disp('unrotated (not to scale)')
% 
% % unit normal
% 
% nhat = [0 1]; % unit normal to average fault plane
% [n(:,1),n(:,2)] = rotate_nt2xy_vec(n(:,1),n(:,2),nhat);
  
% plot profile to scale

% subplot(2,1,2)
% plot(x,y,'b-',xend,yend,'ko')
% % hold on,quiver(x,y,n(:,1)',n(:,2)',0,'r'),hold off
% axis image
% 
% % power spectral density estimators
% 
% figure(2),clf
% clf
% S = spectrum.mtm;
% Hpsd = psd(S,y,'Fs',1/h,'SpectrumType','twosided');
% loglog(k/(2*pi),PSDy,'r',Hpsd.Frequencies,Hpsd.Data,'b',k/(2*pi),PSDy_exact,'k'),xlim([1e-3 max(k/(2*pi))])


% text(0.1,1e-10,'$\beta V_x/\mu$ lalalal','FontName', 'Times', 'FontSize', 25);
% xlabel('$\frac{1}{\alpha}$','FontName', 'Times', 'FontSize', 25);
% set(gca, 'FontName', 'Times', 'FontSize', 25);


