% 1D single sediment layer Amplification Ratio 
% See Papageorgiou and Kim (EESD 1993), Figure 2
% H.Rendon Set. 2001 

home
% Define Model Parameters
a1 = 1775.; b1 = 1025.; ro1 = 1800;
a2 = 4000.; b2 = 2300.; ro2 = 2200;
H = 300;
pi   = acos(-1.);
im   = sqrt(-1.);


% Define the incident wave
theta = input('Enter incidence angle => ');      %Angle of incidence
index = input('(1) P-Wave,    (2) S-Wave   ');
if index == 1           
    p  =  sin(theta*pi/180.)/a2;                %Ray Parameter  P-wave
    fO    = a1/(4.*H);                          %Normalized Frequency for S-wave
else
    p  =  sin(theta*pi/180.)/b2;                %Ray Parameter S-wave
    fO    = b1/(4.*H);                          %Normalized Frequency for S-wave
end
psi2  = sqrt(1./a2^2 - p^2);                    %P-wave vertical slowness Basement
eta2  = sqrt(1./b2^2 - p^2);                    %S-wave vertical slowness Basement
psi1  = sqrt(1./a1^2 - p^2);                    %P-wave vertical slowness Sediment
eta1  = sqrt(1./b1^2 - p^2);                    %S-wave vertical slowness Sediment


fmin = 0.1;       %Minimum frequency
fmax = 10*fO;      %Maximun frequency
N    = 2000;       %Number of frequencies

for i=1:N       %Define counter for the frequencies
  w = 2*pi*( fmin + (i-1)*(fmax-fmin)/(N-1) );
  
  lambda2(1,1) = exp(+im*w*psi2*H);   
  lambda2(2,2) = exp(+im*w*eta2*H);
  lambda2(3,3) = exp(-im*w*psi2*H);
  lambda2(4,4) = exp(-im*w*eta2*H);
  
  e2(1,1)=a2*p;     e2(1,2)=b2*eta2;  e2(1,3)=a2*p;      e2(1,4)=b2*eta2;
  e2(2,1)=a2*psi2;  e2(2,2)=-b2*p;    e2(2,3)=-a2*psi2;  e2(2,4)=b2*p;
  e2(3,1)=2*im*w*ro2*a2*b2^2*p*psi2;  e2(3,2)=im*w*ro2*b2*(1-2*b2^2*p^2);
  e2(3,3)=-2*im*w*ro2*a2*b2^2*p*psi2; e2(3,4)=-im*w*ro2*b2*(1-2*b2^2*p^2);
  e2(4,1)=im*w*ro2*a2*(1-2*b2^2*p^2); e2(4,2)=-2*im*w*ro2*b2^3*p*eta2;
  e2(4,3)=im*w*ro2*a2*(1-2*b2^2*p^2); e2(4,4)=-2*im*w*ro2*b2^3*p*eta2;
  
  e1(1,1)=a1*p;     e1(1,2)=b1*eta1;  e1(1,3)=a1*p;      e1(1,4)=b1*eta1;
  e1(2,1)=a1*psi1;  e1(2,2)=-b1*p;    e1(2,3)=-a1*psi1;  e1(2,4)=b1*p;
  e1(3,1)=2*im*w*ro1*a1*b1^2*p*psi1;  e1(3,2)=im*w*ro1*b1*(1-2*b1^2*p^2);
  e1(3,3)=-2*im*w*ro1*a1*b1^2*p*psi1; e1(3,4)=-im*w*ro1*b1*(1-2*b1^2*p^2);
  e1(4,1)=im*w*ro1*a1*(1-2*b1^2*p^2); e1(4,2)=-2*im*w*ro1*b1^3*p*eta1;
  e1(4,3)=im*w*ro1*a1*(1-2*b1^2*p^2); e1(4,4)=-2*im*w*ro1*b1^3*p*eta1;
  
  lambda1(1,1) = 1.;   
  lambda1(2,2) = 1.;
  lambda1(3,3) = 1.;
  lambda1(4,4) = 1.;
  
  lambda1I(1,1) = exp(-im*w*psi1*H);   
  lambda1I(2,2) = exp(-im*w*eta1*H);
  lambda1I(3,3) = exp(+im*w*psi1*H);
  lambda1I(4,4) = exp(+im*w*eta1*H);
  
  L=e1*lambda1*lambda1I*inv(e1)*e2*lambda2;
    
  A=L(1:2,1:2);     B=L(3:4,1:2);
  if index == 1
      C = .5*L(1:2,3); D=.5*L(3:4,3);
  else
      C = .5*L(1:2,4); D=.5*L(3:4,4);
  end  
  DESP = C - A*inv(B)*D;
  UX(i) = abs(DESP(1)); UZ(i) = abs(DESP(2));
  f(i) = w/2/pi/fO;
  
end

plot(f,UX,f,UZ)
legend('UX','UZ')
xlabel('Normalized Frequency,  f/fO')
ylabel('Amplitude')