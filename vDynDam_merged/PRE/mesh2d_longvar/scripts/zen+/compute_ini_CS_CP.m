% Script to compute initial strain based on fault angle with principal
% stresses as well as the others necessary input parameters

% Marion Thomas, last modified February 2018

% Script to compute initial strain based on fault angle with principal
% stresses as well as the others necessary input parameters

% Marion Thomas, last modified June 2017

%CALLS: 

%==========================================================================
%% A,B,C,O PARAMETERS
clc

D=0.000001;

%Compute the KI parameters needed
[A,B,C,O]= KI_parameters(D,D,alpha,fs);

Dvar = sqrt((pi*D0(1)*(1-nu))/(alpha^3));

%% cs cp FOR REGIME 2

%Variables for Regime 2
A1 = A*Dvar;
B1 = B*Dvar;
Gamma = (3*(1-2*nu))/(2*(1+nu)) + (3*B1^2*(1-2*nu))/(4*(1+nu)) + (A1^2)/2;
% a1 = (1/Gamma)*(1+(B1^2)/2);
% b1 = (-1/Gamma)*((A1*B1)/2);
% b2 = (1/Gamma)*((A1^2)/2 + (3*(1-2*nu))/(2*(1+nu)));

%equivalent Lame parameter
mus     = (1.0/2.0)*(mu/Gamma)*(3.0*(1.0-2.0*nu)/(1+nu) + A1^2);
lambdas = (mu/Gamma)*((3.0*nu)/(1.0+nu) + (B1^2.0)/2.0 - (A1^2)/3.0);
            
%apparent S and P waves velocity
css2 = sqrt(mus/rho);
cps2 = sqrt((lambdas+2*mus)/rho);

%reduction
csred2=100*(cs-css2)/cs
cpred2=100*(cp-cps2)/cp

%% cs cp FOR REGIME 3



%Variables for Regime 3
C1 = C*Dvar;
O1 = O*Dvar;
% C1b2 = ((3*(1-2*nu))/(1+nu)+C1^2)^(-1);
% O1b2 = (2+O1^2)^(-1);
cst1 = ((3.0*(1-2*nu))/(1+nu)+C1^2)^(-1);
cst2 = (2.0+O1^2)^(-1);

%equivalent Lame parameter
mus     = 2.0*mu*cst2;
lambdas = mu*(2.0*cst1-4.0*cst2/3.0);
            
%apparent S and P waves velocity
css3 = sqrt(mus/rho);
cps3 = sqrt((lambdas+2*mus)/rho);
        
%reduction
csred3=100*(cs-css3)/cs
cpred3=100*(cp-cps3)/cp
        
        
        
        
        
        
        
        
