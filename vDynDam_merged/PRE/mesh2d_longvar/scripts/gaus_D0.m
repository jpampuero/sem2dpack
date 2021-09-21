% Script to compute initial damage density (D0) with a gaussian
% distribution around the fault

% Marion Thomas, last modified June 2018

%CALLS: 

%==========================================================================
%% DYNAMIC DAMAGE PARAMETERS

if gaus_tag==1
    
    %end of the exponential profil for the TOP
    if strcmp(distD0T,'pro_zone')==1
        compute_ini_param; clc
        zeT=pro_zone;
    else
        zeT=str2double(distD0T);
    end
    
    %end of the exponential profil for the BOTTOM
    if strcmp(distD0B,'pro_zone')==1
        compute_ini_param; clc
        zeB=pro_zone;
    else
        zeB=str2double(distD0B);
    end
    
    %Set the default values
    D0p=[];D0pz=[];
    zi=0; 
    dz=round(150/h)*h;         %steps between differentes values of D0
    if dz/h <2, dz=2*h; end
    

    %TOP part
    if zeT==zi
        D0expT=[];
        zT=[];
    else
        zT=zi:dz:zeT;%+dz;
        y=exp(-zT.^2/(zeT^2));
        my=min(y);
        D0expT=[y-my]./max(y-my).*(D0iniT-D0endT)+(D0endT);
    end
    
    %BOTTOM part
    if zeB==zi
        D0expB=[];
        zB=[];
    else
        zB=zi:dz:zeB;%+dz;
        y=exp(-zB.^2/(zeB^2));
        my=min(y);
        D0expB=[y-my]./max(y-my).*(D0iniB-D0endB)+(D0endB);
    end

    %Combine values
    D0p=[flip(D0expB(1:end-1))';D0expT(1:end-1)'];
    D0pz=[flip(zB(2:end)')*-1;zT(2:end)'];
%     
%     
% 
% D0p=[0.35;0.3;0.25];D0pz=[0.6e3;1.2e3;1.8e3];

    
end