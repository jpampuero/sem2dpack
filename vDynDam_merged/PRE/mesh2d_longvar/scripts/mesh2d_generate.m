% ========================================================================
% Script to build the grid for a single fault in the middle of a uniform
% space with dynamic damage rheology
% with this script, you can only subdivide the domaine along the z-axis

% Marion Thomas, last modified August 2018

% CALLS: mesh2d_tfi, mesh2d_plot, mesh2d_merge, mesh2d_write
%=========================================================================
%% DOMAIN GEOMETRY

domainlength = abs(BOX(2)-BOX(1));  % lentgh of the domain
domainheight = abs(BOX(4)-BOX(3));  % height of the domain (perpendicular to the fault)

%% SUB-DOMAINS COORDINATES 

%Remove points if the width of the damage zone is 0
if numel(D0)>1
    idD0=find(D0z);
    D0=D0(idD0);D0z=D0z(idD0);
end

% x coordinates of the points on the boundaries of the domain
pxB = BOX(1:2);

% z coordinates of the points on the boundaries of the domain
pzB=[sort([D0z(1);D0z(2:end-1)+pFl(2);pFl(2);D0z(end)]),...
    sort([D0z(1);D0z(2:end-1)+pFr(2);pFr(2);D0z(end)])];

%IF the fault start after the left side of the domain
if pFl(1)<Fl(1)
    px=sort([pxB,Fl(1)]);
    pz=[pzB(:,1),sort([D0z(1);D0z(2:end-1)+Fl(2);Fl(2);D0z(end)]),pzB(:,2)];
else
    px=pxB;pz=pzB;
end

%IF the fault stop before the right side of the domain
if pFr(1)>Fr(1)
    px=sort([px,Fr(1)]);
    pz=[pz(:,1:end-1),sort([D0z(1);D0z(2:end-1)+Fr(2);Fr(2);D0z(end)]),pz(:,end)];
end

%number of domain
if numel(D0)==1, numDz=numel(D0)+1; else numDz=numel(D0); end
numDx=(numel(px)-1);
numD=numDx*numDz;

%Resolution correction

%% NUMBERS OF ELEMENTS IN EACH DOMAINS

% number of elements along x : such that horizontal element size = h
NELX=zeros(numDz,numDx);
for i=1:numDx
    NELX(:,i) =abs(px(i+1)-px(i))/h;
end

% number of elements along z :
NELZ=zeros(numDz,numDx);
indNELZ=find([pFl(2),pFr(2)]==max([pFl(2),pFr(2)]),1);
for i=1:numDz
    NELZ(i,:) =abs(pz(i+1,indNELZ)-pz(i,indNELZ))/h;
end

%Initial damage profil
check_D0

%% SUB-DOMAINS MESHING

pdx=zeros(4,numD);
pdz=zeros(4,numD);

for j=1:numDx
    
    for i=1:numDz
        
        %Sub-domain number
        k=(j-1)*numDz+i;
    
        %x-coordinates of the subdomain
        pdx(1,k) = px(j);	%bottom left
        pdx(2,k) = px(j+1);	%bottom right
        pdx(3,k) = px(j+1);	%top right
        pdx(4,k) = px(j);	%top left
        
        %z-coordinates of the subdomain
        pdz(1,k) = pz(i,j);     %bottom left
        pdz(2,k) = pz(i,j+1);	%bottom right
        pdz(3,k) = pz(i+1,j+1);	%top right
        pdz(4,k) = pz(i+1,j);	%top left
        
        %sides
        dom(k).curves{2} = sample_segments([[pdx(2,k),pdz(2,k)];[pdx(3,k),pdz(3,k)]],NELZ(i,j)); % right 
        dom(k).curves{4} = sample_segments([[pdx(1,k),pdz(1,k)];[pdx(4,k),pdz(4,k)]],NELZ(i,j)); % left  

        %-- mesh the bottom
        if i == 1
        dom(k).curves{1} = sample_segments([[pdx(1,k),pdz(1,k)];[pdx(2,k),pdz(2,k)]],NELX(i,j));   % bottom 
        idx1=find(x==pdx(4,k));
        idx2=find(x==pdx(3,k));
        dz=z(idx1)-pdz(4,k);
        dom(k).curves{3} = [x(idx1:idx2)' z(idx1:idx2)'-dz]; %top                                                     % top
            
        %-- mesh the top
        elseif i == numDz
        idx1=find(x==pdx(1,k));
        idx2=find(x==pdx(2,k));
        dz=z(idx1)-pdz(1,k);
        dom(k).curves{1} = [x(idx1:idx2)' z(idx1:idx2)'-dz]; %bottom                                                       % bottom 
        dom(k).curves{3} = sample_segments([[pdx(4,k),pdz(4,k)];[pdx(3,k),pdz(3,k)]],NELX(i,j));   % top

        %-- mesh the middle
        elseif k>1 && k< numD
        idx1=find(x==pdx(1,k));
        idx2=find(x==pdx(2,k));
        dz=z(idx1)-pdz(1,k);
        dom(k).curves{1} = [x(idx1:idx2)' z(idx1:idx2)'-dz]; %bottom                                                      % bottom 
        idx3=find(x==pdx(4,k));
        idx4=find(x==pdx(3,k));
        dz=z(idx3)-pdz(4,k);
        dom(k).curves{3} = [x(idx3:idx4)' z(idx3:idx4)'-dz]; %top                                                     % top

        else 
            disp('problem with the meshing')
        end
    
    % Generate the mesh
    domain(k) = mesh2d_tfi(dom(k).curves,4);

    % tag the domain
%     if numel(D0)==1 & numel(rho)==1,domain(k).etag(:) = 1;else, domain(k).etag(:) = i;end
    if numel(D0)==1,domain(k).etag(:) = 1;else, domain(k).etag(:) = i;end
    
    end

end

%% PLOT THE SUB-DOMAINS

%plot
T_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape'); clf;
hold on; axis equal

%domain
subplot(2,1,1); hold on
mesh2d_plot_vmar(domain,2)
ylim([BOX(3) BOX(4)])
xlim([BOX(1) BOX(2)])
title('sub-domains')


%% MERGE THE DOMAINS INTO A SINGLE MESH

% Define how to connect the boundaries of the different domains
mtab=zeros(4,numD);

for j=1:numDx
    
    for i=1:numDz
        
        %Sub-domain number
        k=(j-1)*numDz+i;
    
        %left boundary of a subdomain
        bdcoor=dom(k).curves{4};
        if j==1, mtab(4,k)=-4; plot(bdcoor(:,1),bdcoor(:,2),'c-','LineWidth',3);...
        else, mtab(4,k)=k-numDz; plot(bdcoor(:,1),bdcoor(:,2),'k-','LineWidth',3);end

        %right boundary of a subdomain
        bdcoor=dom(k).curves{2};
        if j==numDx, mtab(2,k)=-2; plot(bdcoor(:,1),bdcoor(:,2),'r-','LineWidth',3);...
        else, mtab(2,k)=k+numDz; plot(bdcoor(:,1),bdcoor(:,2),'k-','LineWidth',3);end
        
        %Top boundary of a subdomain
        bdcoor=dom(k).curves{3};
        if k==j*numDz; mtab(3,k)=-3; plot(bdcoor(:,1),bdcoor(:,2),'g-','LineWidth',3);...
        else, mtab(3,k)=k+1; plot(bdcoor(:,1),bdcoor(:,2),'k-','LineWidth',3);end

        %bottom boundary of a subdomain
        bdcoor=dom(k).curves{1};
        if k==(j-1)*numDz+1; mtab(1,k)=-1;plot(bdcoor(:,1),bdcoor(:,2),'b-','LineWidth',3);...
        else, mtab(1,k)=k-1; plot(bdcoor(:,1),bdcoor(:,2),'k-','LineWidth',3);end
        
        %does the bottom boundary corresponds to the top of the fault
        if ((pdx(2,k)-Fr(1)+pdz(2,k)-Fr(2)+...
                pdx(1,k)-Fl(1)+pdz(1,k)-Fl(2))==0)
            mtab(1,k)=-6;
        end

        %does the top boundary corresponds to the bottom of the fault
        if ((pdx(3,k)-Fr(1)+pdz(3,k)-Fr(2)+...
                pdx(4,k)-Fl(1)+pdz(4,k)-Fl(2))==0)
            mtab(3,k)=-5;
        end        

    end
        
end

% merge
mesh = mesh2d_merge(domain,mtab);

%% PLOT THE MESH

%final mesh
subplot(2,1,2); hold on;
mesh2d_plot_vmar(mesh,2)
title('final mesh')
ylim([BOX(3) BOX(4)])
xlim([BOX(1) BOX(2)])

%nucleation prone patch
plot([nuclocX1 nuclocX nuclocX2],[nuclocZ1 nuclocZ nuclocZ2],'-co','markersize',10,'MarkerFaceColor','c','linewidth',2)
