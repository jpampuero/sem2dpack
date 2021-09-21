% Script to compute initial damage density (D0) with a expontial decrease
% away from the fault

% Marion Thomas, last modified June 2018

%CALLS: 

%==========================================================================
%%

pattern = 'pro_zone';

if exp_tag==1
    
    %end of the exponential profil for the TOP
    if endsWith(distD0T,'pro_zone')==1
        compute_ini_param; clc
        zeT=eval(distD0T);
        if Nmat==2, zeT=zeT(2);end
    else
        zeT=str2double(distD0T);
    end
    
    %end of the exponential profil for the BOTTOM
    if endsWith(distD0B,'pro_zone')==1
        compute_ini_param; clc
        zeB=eval(distD0B);
        if Nmat==2, zeB=zeB(1);end
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
        y=exp(-zT/zeT);
        my=min(y);
        D0expT=[y-my]./max(y-my).*(D0iniT-D0endT)+(D0endT);
    end
    
    %BOTTOM part
    if zeB==zi
        D0expB=[];
        zB=[];
    else
        zB=zi:dz:zeB;%+dz;
        y=exp(-zB/zeB);
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