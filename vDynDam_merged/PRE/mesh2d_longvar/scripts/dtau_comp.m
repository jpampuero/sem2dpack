% Script to define de dtau for the nucleation prone patch

% Marion Thomas, Harsha Bhat, last modified June 2018

%CALLS: 

%==========================================================================
%% DEFINE THE DTAU FOR THE NUCLEATION PRONE PATCH

%compute traction
nx = n(:,1);
nz = n(:,2);
Tx = Sxx*nx + Sxz*nz;
Tz = Szz*nz + Sxz*nx;
Tn = Tx.*nx + Tz.*nz;
Tt = Tx.*nz - Tz.*nx;

%nucleation-prone patch
[~,ileft] = find(xF>=nuclocX-nucsize,1);
[~,iright] = find(xF>=nuclocX+nucsize,1);

%Find the value to be added in par.inp for the nucleation patch
fstrength=(MuS*correctionfactor*(-Tn(ileft:iright))-Tt(ileft:iright));
dTt_nuc = max(fstrength);
Tt_nuc = Tt(ileft:iright)+dTt_nuc*ones(length(xF(ileft:iright)),1);
disp(' '); disp(['dtau for the nucleation prone patch : ', num2str(dTt_nuc/1e6),' MPa'])


%% PLOT

if exist('rms')==1
    
    %Plot the fault with tractions vectors
    T_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape'); clf;
    subplot(3,1,1)
    plot(xF,zF,'b-')
    hold on,quiver(xF,zF,n(:,1)',n(:,2)',0.15,'r'),hold off
    ylim([abs(max(zF)-min(zF))/2+min(zF)-500 abs(max(zF)-min(zF))/2+min(zF)+500])
    % ylim(BOX(3:4));
    daspect([1 0.5 1])

    %plot and give the change in shear stress needed
    subplot(3,1,2)
    hold on
    plot(xF,Tt,'-+')
    xlabel('$x$','Fontsize',20)
    ylabel('$T_t$','Fontsize',20)
    plot(xF(ileft:iright),Tt(ileft:iright)+dTt_nuc*ones(length(xF(ileft:iright)),1),'-r');
    text(-nucsize/2,0.7*(max(Tt)+dTt_nuc),['$\Delta T_{t,nuc}$ = ' num2str(dTt_nuc,'%4.3e')],'fontsize',18)
    grid on

    %Apparent friction
    subplot(3,1,3)
    plot(xF,Tt./-Tn);
    hold on
    plot(xF(ileft:iright),Tt_nuc./-Tn(ileft:iright),'-r');
    grid on
    xlabel('$x$','Fontsize',20)
    ylabel('$f_0 = T_{t0}/(-T_{n0})$','Fontsize',20)

end