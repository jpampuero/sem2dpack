% rock amplification factor, P-SV, arbitrary incidence
% INPUT: theta = angle respect to X
%        index = incidence (1) P or (2) SV
%        a2,b2 = P and S wave velocities of the rock
% Reference: Aki & Richards, page 140

% Test: figure 5.6 page 142  
% WARNING: their figure is cut before p = 0.2 (at least for the PP)
index = 1 ; a2 = 5000. ; b2 = 3000.;

N = 100 ;
for itheta = 1:N
    
  theta = 0. + 90.*itheta/N ;
  
% convert incidence angle 
% respect to X ---> respect to DOWN
%theta = theta - 90;

if index == 1           
    p  =  sin(theta*pi/180.)/a2;                %Ray Parameter  P-wave
else
    p  =  sin(theta*pi/180.)/b2;                %Ray Parameter S-wave
end
psi2  = sqrt(1./a2^2 - p^2);                    %P-wave vertical slowness Basement
eta2  = sqrt(1./b2^2 - p^2);                    %S-wave vertical slowness Basement

det = 1./( (1./b2^2-2.*p^2)^2 +4.*p^2 * psi2 * eta2 );

if index == 1 
    PP(itheta) = ( -(eta2^2-p^2)^2 +4.*p^2 * psi2 * eta2 )*det;
    PS(itheta)  = ( 4. *a2/b2 *p *psi2 *(eta2^2-p^2) )*det;
else
    SP(itheta)  = ( 4.*b2/a2 *p *eta2 *(1./b2^2-2.*p^2) )*det;
    SS(itheta)  = ((1./b2^2 - 2.*p^2)^2 -4.*p^2 * psi2 * eta2 )*det;
end

pp(itheta) = p*1000. ;

end

if index == 1 
    plot(pp,PP,'o-',pp,PS,'x-')
    legend('PP','PS')
else
    plot(pp,SP,'o--',pp,SS,'x--')
    legend('SP','SS')
end

xlabel('Slowness p')
ylabel('Coefficient')
