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
idD0=find(D0z);
D0=D0(idD0);D0z=D0z(idD0);

%add the fault as a boundary
D0z=sort([D0z;0]);

% x coordinates of the points on the boundaries
% px = [min(x),max(x)];
px = [BOX(1),BOX(2)];

% z coordinates of the points on the boundaries
pz=[[D0z(1);D0z(2:end-1)+pFl(2);D0z(end)],[D0z(1);D0z(2:end-1)+pFr(2);D0z(end)]];

%number of domain
numD=numel(D0);

% number of elements along x : such that horizontal element size = h
NELX=zeros(numD,1);
for i=1:numD
    NELX(i) =ceil(abs(px(2)-px(1))/h);
end

% number of elements along z :
NELZ=zeros(size(D0));
for i=1:numD
    NELZ(i) =round(abs(D0z(i+1)-D0z(i))/h);
end

%Initial damage profil
% check_D0

%% SUB-DOMAINS MESHING

pdx=zeros(4,numD);
pdz=zeros(4,numD);

for i=1:numD
    
        %x-coordinates of the subdomain
        pdx(1,i) = px(1); 	%bottom left
        pdx(2,i) = px(2);   %bottom right
        pdx(3,i) = px(2);   %top right
        pdx(4,i) = px(1);   %top left
        
        %z-coordinates of the subdomain
        pdz(1,i) = pz(i,1);     %bottom left
        pdz(2,i) = pz(i,2);     %bottom right
        pdz(3,i) = pz(i+1,2);   %top right
        pdz(4,i) = pz(i+1,1);	%top left
        
        %sides
        curves{2} = sample_segments([[pdx(2,i),pdz(2,i)];[pdx(3,i),pdz(3,i)]],NELZ(i)); % right 
        curves{4} = sample_segments([[pdx(1,i),pdz(1,i)];[pdx(4,i),pdz(4,i)]],NELZ(i)); % left  

        %-- mesh the middle
        if i>1 && i< numD
        curves{1} = [x' z'+ D0z(i)];                                                    % bottom 
        curves{3} = [x' z'+ D0z(i+1)];                                                  % top

        %-- mesh the bottom
        elseif i==1
        curves{1} = sample_segments([[pdx(1,i),pdz(1,i)];[pdx(2,i),pdz(2,i)]],NELX(i)); % bottom 
        curves{3} = [x' z'+ D0z(i+1)];                                                  % top
        
        %-- mesh the top
        else
        curves{1} = [x' z'+ D0z(i)];                                                    % bottom 
        curves{3} = sample_segments([[pdx(4,i),pdz(4,i)];[pdx(3,i),pdz(3,i)]],NELX(i)); % top

        end
    
    % Generate the mesh
    domain(i) = mesh2d_tfi(curves,4);

    % tag the domain
    domain(i).etag(:) = i;
    
end

%% MERGE THE DOMAINS INTO A SINGLE MESH

% Define how to connect the boundaries of the different domains
mtab=zeros(4,numD);

for j=1:numD
    
    %right and left boundaries of a subdomain always correspond to -2 and -4
	mtab(2,j)=-2; %right
	mtab(4,j)=-4; %left

    %does the bottom boundary corresponds to the top of the fault
    if ((pdx(2,j)-pFr(1)+pdz(2,j)-pFr(2)+...
            pdx(1,j)-pFl(1)+pdz(1,j)-pFl(2))==0)
        mtab(1,j)=-6;
    else
        mtab(1,j)=j-1;
    end

	%does the top boundary corresponds to the bottom of the fault
    if ((pdx(3,j)-pFr(1)+pdz(3,j)-pFr(2)+...
            pdx(4,j)-pFl(1)+pdz(4,j)-pFl(2))==0)
        mtab(3,j)=-5;
    else
        mtab(3,j)=j+1;
    end

	%1st domain
    if j==1;mtab(1,j)=-1;end
    
    %last domain
    if j==numD; mtab(3,j)=-3;end
        
end

% merge
mesh = mesh2d_merge(domain,mtabx);

%% PLOT THE MESH

%plot
T_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape'); clf;
hold on; axis equal
mesh2d_plot_vmar(domain,2)
% mesh2d_plot(domain,2)
skip = 1;
% plot(xF(1:skip:end),zF(1:skip:end),'r','linewidth',3)
ylim([BOX(3) BOX(4)])
xlim([BOX(1) BOX(2)])