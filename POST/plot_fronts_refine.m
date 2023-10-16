% PLOT_FRONTS space-time plot of rupture front and process zone tail
%
% SYNTAX	[Trup,Tpz] = plot_fronts(V,Veps,D,Dc,X,DT)
%
% INPUTS	V(:,:)	space-time slip velocity data
%		Veps	slip velocity threshold to define rupture front
%		D(:,:)	space-time slip data
%		Dc	slip threshold to define the tail of process zone
%		X(:)	positions of slip and velocity data
%		DT	timestep of slip and velocity data
%		
% OUTPUTS	Trup(:) rupture times (for which V=Veps) at each point 
%			with quadratic time interpolation
%		Tpz(:)  time of tail of process zone (for which D=Dc) 
%			at each location, with quadratic time interpolation
%		Plot  	a shaded area in space-time indicating the head and 
%			tail of the rupture front
%
function [Trup,Tpz] = plot_fronts_refine(V,Veps,D,Dc,X,DT)

%		T(:)	times of slip and velocity data

NX = length(X);
[N1,N2] = size(V);
if N1==NX
  NT=N2;
elseif N2==NX
  V = V';
  D = D';
  NT=N1;
else
  error('Size mismatch (V/D/X)')
end

Trup = repmat(NT,NX,1);
Tpz  = repmat(NT,NX,1);

% force V(:,1)=0 at t=0
%if nnz(V(:,1))
%  V = [zeros(NX,1) V];
%  D = [zeros(NX,1) D];
%end

for k=1:NX,

 % rupture front
  m = find( V(k,:)>Veps );
  if ~isempty(m)
    m = m(1)-1; 
    if m>1
      Trup(k) = m-1 + (Veps-V(k,m))/(V(k,m+1)-V(k,m)); % assumes m=1 is t=0
    else
      Trup(k) = 0;
    end
  end

 % end of process zone
  m = find( D(k,:)<Dc);
  m = m(end);
  if m<NT
    Tpz(k) = m-1 + (Dc-D(k,m))/(D(k,m+1)-D(k,m)); % assumes m=1 is t=0
  end

end

Trup = Trup*DT;
Tpz = Tpz*DT;

%% subplot 1, slip rate color figure + rupture fronts and tails
% Figure size setting
f = figure;
x0=190;
y0=190;
width=800;
height=350;
set(f,'position',[x0,y0,width,height],'papersize',[42,20]);

% make non-linear cpt
V(V < Veps) = 0; % make V(V<Veps)=0 to increase contrast
cMap = flip(pink(256));
dataMax = max(max(V));
dataMin = min(min(V));
centerPoint = dataMax*1/20 + dataMin*19/20;
scalingIntensity = 6;

x = 1:length(cMap); 
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity * x/max(abs(x));

x = sign(x).* exp(abs(x));
x = x - min(x); x = x*511/max(x)+1; 
newMap = interp1(x, cMap, 1:512);


% plot
subplot(1,2,1);
nt=size(V);

h=pcolor(X/1e3,(0:nt(2)-1)*DT,V');
set(h, 'EdgeColor', 'none')
colormap(newMap);
c = colorbar;
c.Label.FontSize = 17;
c.FontSize = 17;
% c.Position = [0.88 0.2 0.025 0.44];
c.Label.String = 'Slip rate (m/s)';

% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',17);
% set(gca,'XTickLabelMode','auto');
title("Slip rate & process zone");
hold on

% add rupture fronts and tails (process zone)
if nargout==0
  % h1=figure;
  hold on
  plot(X/1e3,Tpz,'blue','linewidth',2);
  plot(X/1e3,Trup,'red','linewidth',2);
  %area( X,Tpz,'FaceColor','b')
  %hold on; area( X,Trup,'FaceColor','w'); hold off
  xlabel('Distance (km)','FontSize',17,'FontName', 'Helvetica')
  ylabel('Time (s)','FontSize',17,'FontName', 'Helvetica')
end

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',17)
set(gca,'XTickLabelMode','auto')


%% subplot 2, Vr v.s. distance along fault
subplot(1,2,2);

% compute the rupture speed, below parts are added by Pablo & Liuwei, Sep 19, 2023

% % Option 1: compute speed by interpolation
% x_in = linspace(min(X), max(X), length(X));
% t_in = interp1(X, Trup, x_in);
% vrup = zeros(size(Trup));
% vrup(1) = (x_in(2)-x_in(1))./(t_in(2)-t_in(1));
% vrup(2:end-1) = (x_in(3:end)-x_in(1:end-2))./(t_in(3:end)-t_in(1:end-2)); % central difference
% vrup(end) = (x_in(end)-x_in(end-1))./(t_in(end)-t_in(end-1));
% % end of option 1

% Option 2: compute speed by original grids
dx = diff(X);
dx = [dx(1); dx(:); dx(end)]; % trick to handle the two ends
dt = diff(Trup);
dt = [dt(1); dt(:); dt(end)];
t_in = Trup;
vrup  = ( dx(1:end-1) + dx(2:end) ) ./ ( dx(2:end)./dx(1:end-1).*dt(1:end-1) + dx(1:end-1)./dx(2:end).*dt(2:end) );
% end of option 2

plot(X/1e3, smoothdata(vrup,'gaussian',20),"LineWidth",2,"Color",'black');
hold on;

% add references speeds (Cs and Cp in Par.inp)
Cs = 3330; % m/s, shear wave speed
Cp = 5770; % m/s, P wave speed

% add references speeds (Cs and Cp in Par.inp)
plot(X/1e3,Cs.*ones(1,length(X)),'m--','linewidth',1.5);
plot(X/1e3,Cp.*ones(1,length(X)),'m--','linewidth',1.5);
text((2.2*min(X)+7.7*max(X))/10e3,Cs*1.04,...
    strcat('Cs = ',string(Cs/1000),' km/s'),'FontSize',19,'Color','m');
text((2.2*min(X)+7.7*max(X))/10e3,Cp*1.04,...
    strcat('Cp = ',string(Cp/1000),' km/s'),'FontSize',19,'Color','m');

% refine the figure
xlabel("Distance (km)",'FontName', 'Helvetica');
ylabel("Rupture Speed (m/s)",'FontName', 'Helvetica')
ylim([0,Cp*1.1]);
xlim([min(X/1e3),max(X/1e3)]);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',17)
set(gca,'XTickLabelMode','auto')
title("Speed v.s. distance");

%% plot a subshear case and a supershear case for references
% you can comment the whole part if you do not need references
fileID = fopen('Vr_dist_sub.txt','r');
formatSpec = '%f %f\n';
sizeVr_sub = [2 Inf];
Vr_sub = fscanf(fileID,formatSpec,sizeVr_sub);
fclose(fileID);

fileID = fopen('Vr_dist_super.txt','r');
formatSpec = '%f %f\n';
sizeVr_super = [2 Inf];
Vr_super = fscanf(fileID,formatSpec,sizeVr_super);
fclose(fileID);

plot(Vr_sub(1,:),smoothdata(Vr_sub(2,:),'gaussian',20),'color',[0.75 0.75 0.75],'linewidth',1.8);
plot(Vr_super(1,:),smoothdata(Vr_super(2,:),'gaussian',20),'color',[0.75 0.75 0.75],'linewidth',1.8);
text(37,2600,['Persistent',newline,' subshear'],'FontSize',18,'Color',[0.4 0.4 0.4]);
text(38,5400,['Persistent',' supershear'],'FontSize',19,'Color',[0.4 0.4 0.4]);
text(75,4000,['Intermittent',newline,'supershear'],'FontSize',19,'Color',[0.0 0.0 0.0]);
%% print out

print('Summary_2.pdf', '-dpdf', '-bestfit');
print('Summary_2.png', '-dpng', '-r300');

%% wrtie a txt file of Vr v.s. dist for plotting
fileID = fopen('Vr_dist.txt','w');
for i=1:length(X)
    fprintf(fileID,'%6.2f %6.2f\n',X(i)/1e3,smoothdata(vrup(i),'gaussian',20));
end
fclose(fileID);