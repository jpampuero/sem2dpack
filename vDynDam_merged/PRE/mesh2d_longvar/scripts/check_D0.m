% Script to check the initial damage profil

% Marion Thomas, last modified June 2018

%CALLS: 

%==========================================================================
%%
subplot(1,4,1); hold on

if numel(D0)==1,D0plot=[D0,D0]; else, D0plot=D0; end

for i=1:numDz
    plot(ones(NELZ(i,indNELZ)+1,1)*D0plot(i),[pz(i,indNELZ):h:pz(i+1,indNELZ)]'./min(pro_zone),'*-')
end

ylim([BOX(3)-dz BOX(4)+dz]./pro_zone);
xlim([0 1]); 
ylabel('Distance away from the fault / R_0')
xlabel('D_0')
text(0.7,(BOX(3))./min(pro_zone),['R_0=',num2str(round(min(pro_zone))),' m'])
title('Initial Damage')