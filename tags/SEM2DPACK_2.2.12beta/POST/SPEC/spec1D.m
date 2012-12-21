% 1D amplification spectra
% H. Rendon and J.P. Ampuero - Set. 2001
% INPUT: depthfile, theta, index, velocity model, fmin,fmax,N
% USES: get_sediment

% Define constants
pi   = acos(-1.);
im   = sqrt(-1.);

%Define the diferent depths along the Valley
get_sediment ;

M = size(thickness,1);
H = -thickness(1:M,2);
LOC = thickness(1:M,1)*xscal ;  %Rescaling abscissa

% convert incidence angle 
% respect to X ---> respect to DOWN
inc = theta - 90;

% get the rock amplification factor, for stations outside sediment:
[UXrock,UZrock] = rock_ampli(inc,index,a2,b2);

if index == 1           
    p  =  sin(inc*pi/180.)/a2;                %Ray Parameter  P-wave
    fO    = a1/(4.*H);                          %Normalized Frequency for P-wave
else
    p  =  sin(inc*pi/180.)/b2;                %Ray Parameter S-wave
    fO    = b1/(4.*H);                          %Normalized Frequency for S-wave
end
psi2  = sqrt(1./a2^2 - p^2);                    %P-wave vertical slowness Basement
eta2  = sqrt(1./b2^2 - p^2);                    %S-wave vertical slowness Basement
psi1  = sqrt(1./a1^2 - p^2);                    %P-wave vertical slowness Sediment
eta1  = sqrt(1./b1^2 - p^2);                    %S-wave vertical slowness Sediment

for i=1:N       %Define counter for the frequencies
    
  f(i) = fmin + (i-1)*(fmax-fmin)/(N-1) ;
  w = 2*pi*f(i) ;
  iw = im*w ;
  
  e2(1,1)=a2*p;     e2(1,2)=b2*eta2;  
  e2(1,3)=e2(1,1);  e2(1,4)=e2(1,2);
  e2(2,1)=a2*psi2;  e2(2,2)=-b2*p;    
  e2(2,3)=-e2(2,1); e2(2,4)=-e2(2,2);
  e2(3,1)=2*iw*ro2*a2*b2^2*p*psi2;  e2(3,2)=iw*ro2*b2*(1-2*b2^2*p^2);
  e2(3,3)=-e2(3,1);                 e2(3,4)=-e2(3,2) ;
  e2(4,1)=iw*ro2*a2*(1-2*b2^2*p^2); e2(4,2)=-2*iw*ro2*b2^3*p*eta2;
  e2(4,3)=e2(4,1) ;                 e2(4,4)=e2(4,2);
  
  e1(1,1)=a1*p;     e1(1,2)=b1*eta1;  
  e1(1,3)=e1(1,1);  e1(1,4)=e1(1,2);
  e1(2,1)=a1*psi1;  e1(2,2)=-b1*p;    
  e1(2,3)=-e1(2,1); e1(2,4)=-e1(2,2);
  e1(3,1)=2*iw*ro1*a1*b1^2*p*psi1;  e1(3,2)=iw*ro1*b1*(1-2*b1^2*p^2);
  e1(3,3)=-e1(3,1);                 e1(3,4)=-e1(3,2) ;
  e1(4,1)=iw*ro1*a1*(1-2*b1^2*p^2); e1(4,2)=-2*iw*ro1*b1^3*p*eta1;
  e1(4,3)=e1(4,1) ;                 e1(4,4)=e1(4,2);
  
  L1=e1;
  L2=inv(e1)*e2 ;
  
for j=1:M       %Define counter for different depths
  
  if H(j) < 0.  % sediment
      
    iwH = iw*H(j) ;

    lambda2(1,1) = exp(+iwH*psi2);   
    lambda2(2,2) = exp(+iwH*eta2);
    lambda2(3,3) = exp(-iwH*psi2);
    lambda2(4,4) = exp(-iwH*eta2);
 
    lambda1I(1,1) = exp(-iwH*psi1);   
    lambda1I(2,2) = exp(-iwH*eta1);
    lambda1I(3,3) = exp(+iwH*psi1);
    lambda1I(4,4) = exp(+iwH*eta1);
  
    L=L1*lambda1I*L2*lambda2;
    
    A=L(1:2,1:2);     B=L(3:4,1:2);
    if index == 1
      C = .5*L(1:2,3); D=.5*L(3:4,3);
    else
      C = .5*L(1:2,4); D=.5*L(3:4,4);
    end  
    DESP = C - A*inv(B)*D;
    UX(i,j) = abs(DESP(1)); UZ(i,j) = abs(DESP(2));
    
    
  else  % rock
      
    UX(i,j) = UXrock ; UZ(i,j) = UZrock ;
    
  end
  
  
  
end


end
