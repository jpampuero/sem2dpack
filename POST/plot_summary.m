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
function [Trup,Tpz] = plot_summary(V,Veps,D,Dc,X,DT)

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

subplot(2,2,1);
nt=size(V);
mesh([0:nt(2)-1]*DT,X/1e3,V)  %added by Liuwei, Oct 3, 2023
colorbar;
view(270,90);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',17);
set(gca,'XTickLabelMode','auto');
title("Slip rate");

subplot(2,2,2);

if nargout==0
  % h1=figure;
  plot(X,Tpz,X,Trup);
  %area( X,Tpz,'FaceColor','b')
  %hold on; area( X,Trup,'FaceColor','w'); hold off
  xlabel('X','FontSize',17,'FontName', 'Helvetica')
  ylabel('T','FontSize',17,'FontName', 'Helvetica')
end

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',17)
set(gca,'XTickLabelMode','auto')

title("Process zone");

% %print out
% set(h1,'Units','Inches');
% pos = get(h1,'Position');
% set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% 
% print(h1,'Process_zone.pdf', '-dpdf', '-r0');
% print(h1,'Process_zone.png', '-dpng', '-r300');
%% compute the rupture speed, below parts are added by Pablo & Liuwei, Sep 19, 2023

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

% plot speed v.s. time
% figure; plot(t_in,vrup);

%or apply some smoothing, for instance:
% h=figure; 
subplot(2,2,3);
plot(t_in, smoothdata(vrup,'gaussian',20),"LineWidth",1,"Color",'black');
hold on;

% add references speeds (Cs and Cp in Par.inp)
Cs = 3330; % m/s, shear wave speed
Cp = 5770; % m/s, P wave speed
plot(t_in,Cs.*ones(1,length(t_in)),'red');
plot(t_in,Cp.*ones(1,length(t_in)),'blue');
text((min(t_in)+max(t_in))/2,Cs*1.1,...
    strcat('Cs = ',string(Cs/1000),' km/s'),'FontSize',15,'Color','red');
text((min(t_in)+max(t_in))/2,Cp*1.05,...
    strcat('Cp = ',string(Cp/1000),' km/s'),'FontSize',15,'Color','blue');

% refine the figure
xlabel("Time (s)",'FontName', 'Helvetica');
ylabel("Rupture Speed (m/s)",'FontName', 'Helvetica')
ylim([0,Cp*1.1]);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',17)
set(gca,'XTickLabelMode','auto')
title("Speed v.s. time");


% plot speed v.s. distance
subplot(2,2,4);
plot(X, smoothdata(vrup,'gaussian',20),"LineWidth",1,"Color",'black');
hold on;

% add references speeds (Cs and Cp in Par.inp)
plot(X,Cs.*ones(1,length(X)),'red');
plot(X,Cp.*ones(1,length(X)),'blue');
text((min(X)+max(X))/2,Cs*1.1,...
    strcat('Cs = ',string(Cs/1000),' km/s'),'FontSize',15,'Color','red');
text((min(X)+max(X))/2,Cp*1.05,...
    strcat('Cp = ',string(Cp/1000),' km/s'),'FontSize',15,'Color','blue');

% refine the figure
xlabel("Distance (m)",'FontName', 'Helvetica');
ylabel("Rupture Speed (m/s)",'FontName', 'Helvetica')
ylim([0,Cp*1.1]);
xlim([min(X),max(X)]);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',17)
set(gca,'XTickLabelMode','auto')
title("Speed v.s. distance");

%print out
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print('Summary.pdf', '-dpdf', '-bestfit');
print('Summary.png', '-dpng', '-r300');