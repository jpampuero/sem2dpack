% MESH2D_EX0	2D mesh for a vertical cross-section across a vertical strike-slip fault
clc;close all;clear all
%-- you can modify the following parameters:
% BOX = [-15e3 15e3 -4.5e3 4.5e3]; 	% domain limits [xmin xmax zmin zmax]
% FAULT_X = BOX(2)-6e3;
% FAULT_BOTTOM = -1.5e3; 	% depth of the bottom tip of the fault
% FAULT_TOP = 1.5e3; 	% depth of the bottom tip of the fault
% h = 150;			% element size
BOX = [-9e3 9e3 -3e3 3e3];   % domain limits in m [xmin xmax zmin zmax]
FAULT_X = BOX(2)-4e3;
FAULT_BOTTOM = -1e3; 	% depth of the bottom tip of the fault
FAULT_TOP = 1e3; 	% depth of the bottom tip of the fault
h = 75;			% element size

%-- you don't need to modify anything below this line

% NELX_L = ceil((FAULT_X-BOX(1))/h);
% NELX_R = ceil((BOX(2)-FAULT_X)/h);
% NELZ_T = ceil((FAULT_TOP-FAULT_BOTTOM)/h);
% NELZ_B = ceil((FAULT_BOTTOM-BOX(3))/h);
% if FAULT_TOP~=0
%     NELZ_T2 = ceil((BOX(4)-FAULT_TOP)/h);
% end 
% NELX_L = (FAULT_X-BOX(1))/h
% NELX_R = (BOX(2)-FAULT_X)/h
% NELZ_T = (FAULT_TOP-FAULT_BOTTOM)/h
% NELZ_B = (FAULT_BOTTOM-BOX(3))/h
% if FAULT_TOP~=0
%     NELZ_T2 = (BOX(4)-FAULT_TOP)/h
% end 

resadj=4*h;
% box_length=NELX_L+NELX_R;
% box_width=NELZ_T2+NELZ_B+NELZ_T;
% BOX=[-ceil(box_length/2/9)*9 ceil(box_length/2/9)*9 -ceil(box_width/2/9)*9 ceil(box_width/2/9)*9]*h;
BOX=floor(BOX./resadj)*resadj
FAULT_X = floor((BOX(2)-4e3)/resadj)*resadj;
FAULT_BOTTOM = floor(-1e3/resadj)*resadj; 	% depth of the bottom tip of the fault
FAULT_TOP = floor(1e3/9/h)*resadj; 	% depth of the bottom tip of the fault



NELX_L = ceil((FAULT_X-BOX(1))/h);
NELX_R = ceil((BOX(2)-FAULT_X)/h);
NELZ_T = ceil((FAULT_TOP-FAULT_BOTTOM)/h);
NELZ_B = ceil((FAULT_BOTTOM-BOX(3))/h);
if FAULT_TOP~=0
    NELZ_T2 = ceil((BOX(4)-FAULT_TOP)/h);
end 

% box vertices
p1 = [BOX(1) BOX(3)];
p2 = [BOX(2) BOX(3)];
p3 = [BOX(2) BOX(4)];
p4 = [BOX(1) BOX(4)];

% points on the fault
p7 = [FAULT_X BOX(4)];
p9 = [FAULT_X FAULT_BOTTOM];

% points on box edges, projecting fault bottom tip
p5 = [FAULT_X BOX(3)];
p6 = [BOX(2) FAULT_BOTTOM];
p8 = [BOX(1) FAULT_BOTTOM];
if FAULT_TOP~=0
    p10=p7; p7=[FAULT_X FAULT_TOP];
    p11=p3;p3=[BOX(2) FAULT_TOP];
    p12=p4;p4=[BOX(1) FAULT_TOP];
end

% bottom left domain
ps = [p1; p5; p9; p8];
X=ps(:,1); Y=ps(:,2);
curves{1} = sample_segments([ [X(1) Y(1)]; [X(2) Y(2)] ],NELX_L);
curves{2} = sample_segments([ [X(2) Y(2)]; [X(3) Y(3)] ],NELZ_B);
curves{3} = sample_segments([ [X(4) Y(4)]; [X(3) Y(3)] ],NELX_L);
curves{4} = sample_segments([ [X(1) Y(1)]; [X(4) Y(4)] ],NELZ_B);
% domain(1)=mesh2d_tfi(curves,4);
domain(1) = mesh2d_quad(ps(:,1),ps(:,2),NELX_L,NELZ_B,4);
domain(1).etag(:)=1;
disp(['domain 1 has ',num2str(NELX_L),'x',num2str(NELZ_B),' cells and then /4 it is ',num2str(NELX_L/4),'x',num2str(NELZ_B/4)])

% bottom right domain
ps = [p5; p2; p6; p9];
X=ps(:,1); Y=ps(:,2);
curves{1} = sample_segments([ [X(1) Y(1)]; [X(2) Y(2)] ],NELX_R);
curves{2} = sample_segments([ [X(2) Y(2)]; [X(3) Y(3)] ],NELZ_B);
curves{3} = sample_segments([ [X(4) Y(4)]; [X(3) Y(3)] ],NELX_R);
curves{4} = sample_segments([ [X(1) Y(1)]; [X(4) Y(4)] ],NELZ_B);
domain(4)=mesh2d_tfi(curves,4);
% domain(2) = mesh2d_quad(ps(:,1),ps(:,2),NELX_R,NELZ_B);
domain(4).etag(:)=1;
disp(['domain 4 has ',num2str(NELX_R),'x',num2str(NELZ_B),' cells and then /4 it is ',num2str(NELX_R/4),'x',num2str(NELZ_B/4)])

% top right domain
ps = [p9; p6; p3; p7];
X=ps(:,1); Y=ps(:,2);
curves{1} = sample_segments([ [X(1) Y(1)]; [X(2) Y(2)] ],NELX_R);
curves{2} = sample_segments([ [X(2) Y(2)]; [X(3) Y(3)] ],NELZ_T);
curves{3} = sample_segments([ [X(4) Y(4)]; [X(3) Y(3)] ],NELX_R);
curves{4} = sample_segments([ [X(1) Y(1)]; [X(4) Y(4)] ],NELZ_T);
domain(5)=mesh2d_tfi(curves,4);
% domain(3) = mesh2d_quad(ps(:,1),ps(:,2),NELX_R,NELZ_T);
domain(5).etag(:)=2;
disp(['domain 5 has ',num2str(NELX_R),'x',num2str(NELZ_T),' cells and then /4 it is ',num2str(NELX_R/4),'x',num2str(NELZ_T/4)])

% top left domain
ps = [p8; p9; p7; p4];
X=ps(:,1); Y=ps(:,2);
curves{1} = sample_segments([ [X(1) Y(1)]; [X(2) Y(2)] ],NELX_L);
curves{2} = sample_segments([ [X(2) Y(2)]; [X(3) Y(3)] ],NELZ_T);
curves{3} = sample_segments([ [X(4) Y(4)]; [X(3) Y(3)] ],NELX_L);
curves{4} = sample_segments([ [X(1) Y(1)]; [X(4) Y(4)] ],NELZ_T);
domain(2)=mesh2d_tfi(curves,4);
% domain(4) = mesh2d_quad(ps(:,1),ps(:,2),NELX_L,NELZ_T);
domain(2).etag(:)=2;
disp(['domain 2 has ',num2str(NELX_L),'x',num2str(NELZ_T),' cells and then /4 it is ',num2str(NELX_L/4),'x',num2str(NELZ_T/4)])

if FAULT_TOP~=0
    ps = [p4; p7; p10; p12];
    X=ps(:,1); Y=ps(:,2);
    curves{1} = sample_segments([ [X(1) Y(1)]; [X(2) Y(2)] ],NELX_L);
    curves{2} = sample_segments([ [X(2) Y(2)]; [X(3) Y(3)] ],NELZ_T2);
    curves{3} = sample_segments([ [X(4) Y(4)]; [X(3) Y(3)] ],NELX_L);
    curves{4} = sample_segments([ [X(1) Y(1)]; [X(4) Y(4)] ],NELZ_T2);
    domain(3)=mesh2d_tfi(curves,4);
    %domain(5) = mesh2d_quad(ps(:,1),ps(:,2),NELX_L,NELZ_T2);
    domain(3).etag(:)=3;
    disp(['domain 1 has ',num2str(NELX_L),'x',num2str(NELZ_T2),' cells and then /4 it is ',num2str(NELX_L/4),'x',num2str(NELZ_T2/4)])

    ps = [p7; p3; p11; p10];
    X=ps(:,1); Y=ps(:,2);
    curves{1} = sample_segments([ [X(1) Y(1)]; [X(2) Y(2)] ],NELX_R);
    curves{2} = sample_segments([ [X(2) Y(2)]; [X(3) Y(3)] ],NELZ_T2);
    curves{3} = sample_segments([ [X(4) Y(4)]; [X(3) Y(3)] ],NELX_R);
    curves{4} = sample_segments([ [X(1) Y(1)]; [X(4) Y(4)] ],NELZ_T2);
    domain(6)=mesh2d_tfi(curves,4);
    %domain(6) = mesh2d_quad(ps(:,1),ps(:,2),NELX_R,NELZ_T2);
    domain(6).etag(:)=3;
    disp(['domain 1 has ',num2str(NELX_R),'x',num2str(NELZ_T2),' cells and then /4 it is ',num2str(NELX_R/4),'x',num2str(NELZ_T2/4)])

end

% mesh2d_plot(domain,1)
% axis equal
% ylim([BOX(3) BOX(4)])
% xlim([BOX(1) BOX(2)])

% % tag the elements close to the fault
% % to create a Kelvin-Voigt viscous layer
% e = sub2ind([NELX_R NELZ_B], NELX_L,NELZ_B);
% domain(1).etag(e) = 3;
% e = sub2ind([NELX_R NELZ_B], 1,NELZ_B);
% domain(2).etag(e) = 3;
% e = sub2ind([NELX_R NELZ_T], ones(NELZ_T,1),[1:NELZ_T]');
% domain(3).etag(e) = 3;
% e = sub2ind([NELX_L NELZ_T], NELX_L*ones(NELZ_T,1),[1:NELZ_T]');
% % domain(4).etag(e) = 3;
% 
% merge domains into a single mesh
if FAULT_TOP~=0
%     mtab = [ [ -1 -1  2  1  4  3 ] ; ...
%              [  2 -2 -2 -5  6 -2 ] ; ...
%              [  4  3  6  5 -3 -3 ] ; ...
%              [ -4  1 -6 -4 -4  5 ] ];
%     mtab = [ [ -1 -1  2 -6  4  3 ] ; ...
%              [  2 -2 -2  3  6 -2 ] ; ...
%              [ -5  3  6  5 -3 -3 ] ; ...
%              [ -4  1  4 -4 -4  5 ] ];
    mtab = [ [ -1 -6  2 -1  4  5 ] ; ...
             [  4  5  6 -2 -2 -2 ] ; ...
             [ -5  3 -3  5  6 -3 ] ; ...
             [ -4 -4 -4  1  2  3 ] ];
else
    mtab = [ [ -1 -1  2  1] ; ...
         [  2 -2 -2 -5] ; ...
         [  4  3 -3 -3] ; ...
         [ -4  1 -6 -4] ];
end
    mesh = mesh2d_merge(domain,mtab);


T_handle =figure('Position',[200 500 1600 750],'PaperOrientation','landscape'); clf;
hold on; axis equal
clf
mesh2d_plot_vmar(mesh,2)
axis equal
ylim([BOX(3) BOX(4)])
xlim([BOX(1) BOX(2)])

% %-- export the mesh file
% mesh2d_write(mesh,'ex0');
