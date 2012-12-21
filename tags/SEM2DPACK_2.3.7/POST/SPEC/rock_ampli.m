% [UX,UZ] = rock_ampli(inc,wav,a,b)
%
% PURPOSE	Rock amplification factor (half-space), P-SV, arbitrary incidence
%
% INPUT		theta	incidence angle (degrees with respect to vertical DOWN)
%        	wav 	incidence wave type: 1=P or 2=SV
%        	a,b 	P and S wave velocities of the rock
%
% Reference: Aki & Richards, page 140

function [UX,UZ] = rock_ampli(inc,wav,a,b)

sini = sin(inc*pi/180.);

if wav == 1           
    p  =  sini/a;                %Ray Parameter  P-wave
else
    p  =  sini/b;                %Ray Parameter S-wave
end
psi  = sqrt(1./a^2 - p^2);                    %P-wave vertical slowness Basement
eta  = sqrt(1./b^2 - p^2);                    %S-wave vertical slowness Basement
cosi = a*psi;
cosj = b*eta;
sinj = b*p;

det = 1./( (eta^2-p^2)^2 +4.*p^2 * psi * eta );

if wav == 1 
    PP = ( -(eta^2-p^2)^2 +4.*p^2 * psi * eta )*det;
    PS = ( 4. *a/b *p *psi *(eta^2-p^2) )*det;
    UX = sini*PP + cosj*PS;
    UZ = cosi*PP - sinj*PS;
else
    SP = ( 4.*b/a *p *eta *(eta^2-p^2) )*det;
    SS = ((eta^2-p^2)^2 -4.*p^2 * psi * eta )*det;
    UX = sini*SP + cosj*SS;
    UZ = cosi*SP - sinj*SS;
end
