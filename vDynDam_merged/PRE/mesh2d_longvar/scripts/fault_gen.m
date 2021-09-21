% Script to create the fault for the mesh

% Marion Thomas, last modified August 2018

%CALLS: 

%==========================================================================
%% FAULT GENERATION

if exist('x') == 0
    
    % Flat fault
    if exist('Fl') == 1
        %Correction based on the resolution
        Fl=resocorrec(Fl,h,Q);
        Fr=resocorrec(Fr,h,Q);
        %compute angle
        angF=atand((Fr(2)-Fl(2))/(Fr(1)-Fl(1)));
        nfault=['flatF_ang',num2str(angF)];
        %define the x and z coordinates for the fault
        xF=Fl(1):h:Fr(1);
        zF=Fl(2)+(xF-Fl(1))*sind(angF);
        n = [ones(size(xF,2),1)*sind(angF),ones(size(xF,2),1)*cosd(angF)];

    %self-similar rough fault
    elseif exist('rms') == 1
        [xF,zF,n,s] = rough_flt_gen(rms,seed,BOX,h);
        angF=0;
        
    %default fault
    else
        disp('on va creer une jolie petite faille plate au milieu de la grille')
        xF=BOX(1):h:BOX(2);
        zF=ones(size(xF))*(BOX(3)+(BOX(4)-BOX(3))/2);
        n = [zeros(size(xF,2),1),ones(size(xF,2),1)];
        nfault='flatF';
    end
end

%% FAULT CORRECTION (if necessary)

% xF=resocorrec(xF,h,Q);
% zF=resocorrec(zF,h,Q);
x=xF;
z=zF;

% plot
T_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape'); clf;
subplot(1,4,2:4); hold on;
plot(xF./min(pro_zone),zF./min(pro_zone),'ko')

%Left border correction
if min(xF) ~= BOX(1)
    if min(xF) > BOX(1)
        xsup1=BOX(1):h:xF(1); x=[xsup1(1:end-1) x]; 
        zsup1=tand(angF).*(xsup1-xF(1)+h)+zF(1);
        z=[zsup1(1:end-1)-(zsup1(end)-zF(1)) z];
    end
    if min(xF) < BOX(1)
        id1=find(xF==BOX(1)); xF=xF(id1:end);x=xF; 
        zF=zF(id1:end);n=n(id1:end,:);z=zF;
    end
end

%Right border correction
if max(xF) ~= BOX(2)
    if max(xF) < BOX(2)
        xsup2=xF(end):h:BOX(2); x=[x xsup2(2:end)];
        zsup2=tand(angF)*(xsup2-xF(end)-h)+zF(end);
        z=[z zsup2(2:end)-(zsup2(1)-zF(end))];
    end
    if max(xF) > BOX(2)
        id2=find(xF==BOX(2)); xF=xF(1:id2); zF=zF(1:id2);
        id3=find(x==BOX(2)); x=x(1:id3); z=z(1:id3);n=n(1:id3,:);
    end
end

% %Boundary coordinates correction
% x=resocorrec(x,h,Q);
% z=resocorrec(z,h,Q);

%Find where the fault meet the borders of the domain
if exist('Fl') == 1
    pFl=[BOX(1) round((Fl(2)-tand(angF)*(Fl(1)-BOX(1)))/h)*h]; 
    pFr=[BOX(2) round((Fr(2)+tand(angF)*(BOX(2)-Fr(1)))/h)*h];
else
	pFl=[BOX(1),0];Fl=[BOX(1),0];
	pFr=[BOX(2),0];Fr=[BOX(2),0];
end

%%  FIND THE LOCATION OF THE NUCLEATION PRONE PATCH

indnucZ=find(min(abs(xF-nuclocX))==abs(xF-nuclocX));
% nuclocX=resocorrec(xF(indnucZ),h,Q);
% nuclocZ=resocorrec(zF(indnucZ),h,Q);
nuclocX=xF(indnucZ);
nuclocZ=zF(indnucZ);
indnucZ1=find(min(abs(xF-nuclocX-nucsize))==abs(xF-nuclocX-nucsize));
nuclocX1=xF(indnucZ1);
nuclocZ1=zF(indnucZ1);
indnucZ2=find(min(abs(xF-nuclocX+nucsize))==abs(xF-nuclocX+nucsize));
nuclocX2=xF(indnucZ2);
nuclocZ2=zF(indnucZ2);

%% PLOT TO CHECK

plot(xF./min(pro_zone),zF./min(pro_zone),'*')
plot(x./min(pro_zone),z./min(pro_zone),'r-','linewidth',1.5)
plot([BOX(1) BOX(1:2) flip(BOX(1:2))]./min(pro_zone),...
    [BOX(3:4) flip(BOX(3:4)) BOX(3)]./min(pro_zone),'g-')
plot([nuclocX1 nuclocX nuclocX2]./min(pro_zone),[nuclocZ1 nuclocZ nuclocZ2]./min(pro_zone),'-co','markersize',10,'MarkerFaceColor','c','linewidth',2)
% plot([nuclocX1 nuclocX2]./min(pro_zone),[nuclocZ1 nuclocZ2]./min(pro_zone),'co')

dx=(BOX(2)-BOX(1))/10;
dz=(BOX(4)-BOX(3))/10;

%Legend
legend('initial fault','corrected fault','subdomains border','domain limits','nucleation patch')
title('Domain geometry')
axis equal
xlim([BOX(1)-dx BOX(2)+dx]./min(pro_zone));
ylim([BOX(3)-dz BOX(4)+dz]./min(pro_zone));
xlabel('Distance/R_0')
ylabel('Distance/R_0')

