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
        x=BOX(1):h:BOX(2);
        
        %Left border correction
        if Fl(1) < BOX(1); Fl(1)= BOX(1);'coucou'; end
        zBl=-(Fl(1)-x(1))*tand(angF)+Fl(2);
        z=-(Fl(1)-x)*tand(angF)+Fl(2);

        %Right border correction
        if max(Fr(1))> BOX(2); Fr(1)= BOX(2); Fr(2)=z(end); end
        n = [ones(size(x,2),1)*sind(angF),ones(size(x,2),1)*cosd(angF)];
        
        %Point on the fault
        id1=find(abs(x-Fl(1))==min(abs(x-Fl(1))));
        id2=find(abs(x-Fr(1))==min(abs(x-Fr(1))));
        xF=x(id1:id2); zF=z(id1:id2);


    %self-similar rough fault
    elseif exist('rms') == 1
        [xF,zF,n,s] = rough_flt_gen(rms,seed,BOX,h);
        angF=0;
        
    %default fault
    else
        disp('on va creer une jolie petite faille plate au milieu de la grille')
        xF=BOX(1):h:BOX(2);
        zF=ones(size(xF))*(BOX(3)+(BOX(4)-BOX(3))/2);
        x=xF;z=zF;
        Fl=[x(1) z(1)]; Fr=[x(end) z(end)];
        n = [zeros(size(xF,2),1),ones(size(xF,2),1)];
        nfault='flatF';
    end
end

%Define where the fault meet the borders of the domain
pFl=[x(1) z(1)]; 
pFr=[x(end) z(end)];


%% PLOT TO CHECK

% plot
T_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape'); clf;
subplot(1,4,2:4); hold on;
plot(xF./min(pro_zone),zF./min(pro_zone),'ko')



%plot
plot(xF./min(pro_zone),zF./min(pro_zone),'*')
plot(x./min(pro_zone),z./min(pro_zone),'r-')
plot([BOX(1) BOX(1:2) flip(BOX(1:2))]./min(pro_zone),...
    [BOX(3:4) flip(BOX(3:4)) BOX(3)]./min(pro_zone),'g-')

dx=(BOX(2)-BOX(1))/10;
dz=(BOX(4)-BOX(3))/10;

%Legend
legend('initial fault','corrected fault','subdomains border','domain limits')
title('Domain geometry')
axis equal
xlim([BOX(1)-dx BOX(2)+dx]./min(pro_zone));
ylim([BOX(3)-dz BOX(4)+dz]./min(pro_zone));
xlabel('Distance/R_0')
ylabel('Distance/R_0')

