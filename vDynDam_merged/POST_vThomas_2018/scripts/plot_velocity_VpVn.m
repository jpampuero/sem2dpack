% ************************************************************************
% This function plot the particules velocity for the different "seismometers"
% (point outputs, definied in the input file)

% Marion Thomas, last modified october 2019

% CALLS: export_fig.m;

%% ************************************************************************
function [] = plot_velocity_VpVn(Pcoord,data,fault,tend,save_tag,namef,opt)

%% INPUTS VARIABLES

%Pcoord:	seismographs locations
%data:      seismograms (potentially including several simulations)
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef      field to be plot
%field      definie which component (z or x) should be plotted

%options
proZ = opt(2);	%process zone size
fs = opt(1);	%font size
dur=opt(3);     %max duration to plot the signal

%% PLOT PARAMETERS I

%Launch figure
v_handle =figure('Position',[0 500 1600 1000],'PaperOrientation','landscape');  
fv=v_handle;
set(fv,'PaperOrientation','landscape');
set(fv,'PaperUnits','normalized');
set(fv,'PaperPosition', [0 0 1 1]);
namefile=[namef,'velocity_vpvn'];

%selected seismograms
nsS=size(Pcoord,1);         %number of selected seismograms
uxc=unique(Pcoord(:,1));    %unique xcoord
cf_uc = 2;
co=numel(uxc)*cf_uc;        %number of column
li=ceil(nsS/numel(uxc));    %number of lign
if(mod(co,1) == 0), else, disp('even number, cannot make the figure');exit; end

%curve color
% ctabvn=[1 0 0;0.4 0 0;0.8 0 0;0.2 0 0;0.6 0 0];
% ctabvp=[0 0 1;0 0 0.4;0 0 0.8;0 0 0.2;0 0 0.6];
ctabvn=[1 0 0;0.4 0 0;0 0 0];
ctabvp=[0.3 0.6 1;0 0 1;0 0 0];

%% DOWNLOAD DATA

%number of simulations
Nsimu=size(data,1);
s(1).nt=0;
s(1).ns=0;

for i=1:Nsimu
    
%field to plot
    s(i).vp=data(i).ux;
    s(i).vn=data(i).uz;

%corresponding time
    s(i+1).nt=data(i).nt;       %number of time steps
    s(i).t=fault(i).t;          %associated time
    s(i).tend=tend(i);          %end of simulation

%x and z coordinates
    s(i+1).ns=data(i).nsta;     %number of seismograms
    s(i).x=data(i).x;           %x-coord of seismograms
    s(i).z=data(i).z;           %z-coord of seismograms
    
%fault coordinates    
    f(i).x=fault(i).x;           %x-coord of fault
    f(i).z=fault(i).z;           %z-coord of fault

%% SELECTED SEISMOGRAMS

    %find the unique x coordinates
    xsens=Pcoord(:,1);
    xloc1=unique(Pcoord(:,1));
    %only keep the ones that correspond to a x coordinates one the fault
    xloc=xloc1(find(xloc1>=f(i).x(1) & xloc1 <= f(i).x(end)));
    %find the corresponding z coordinates on the fault
    z_xloc=zeros(size(xloc));
    for j=1:numel(xloc)
        dxloc=abs(xloc(j)-f(i).x);
        z_xloc(j)=f(i).z(find(dxloc==min(dxloc)));
    end
    %Define the z coordinates of the sismometers
    zsens=Pcoord(:,2);
    for j=1:numel(xloc)
        idS=find(xsens == xloc(j));
        zsens(idS)=zsens(idS)+z_xloc(j);
    end

    %Find the seismograms the closest to the location given by the user
    taperdur=s(i).tend-s(i).t(1);
    for k=1:nsS
        matfind_p=abs(s(i).x-xsens(k))+abs(s(i).z-zsens(k));
        s(i).ind(k)= find(min(matfind_p) == matfind_p);  %n0 of the selected seismograms    %min/max

    %min/max
        du=(s(i).vn(2:end,s(i).ind(k))-s(i).vn(1:end-1,s(i).ind(k)));
        du=du/max(du);
        imX=find(du>=0.25,1);
        s(i).mX(k)=round(100*(s(i).t(imX)-taperdur*0.1))/100;
        if s(i).mX(k)+dur <= s(i).tend, s(i).MX(k)=s(i).mX(k)+dur; else...
                s(i).MX(k)= s(i).tend; s(i).mX(k)=s(i).MX(k)-dur; end
        s(i).mY(k)=floor(2*(min(min(min(s(i).vn(:,s(i).ind(k)))),min(min(s(i).vp(:,s(i).ind(k)))))))/2;
        s(i).MY(k)=ceil(2*(max(max(max(s(i).vn(:,s(i).ind(k)))),max(max(s(i).vp(:,s(i).ind(k)))))))/2;
   end
    
end

%% PLOT PARAMETERS II

%declare variable
mY=zeros(1,nsS);
MY=zeros(1,nsS);
mX=zeros(1,nsS);
MX=zeros(1,nsS);

%min/max for amplitude and time
for k=1:nsS
    if size(opt,2) == 7
        mY(k) = opt(4);          %min value for the time
        MY(k) = opt(5);          %max value for the time
        mX(k) = opt(6);          %min value for the amplitude
        MX(k) = opt(7);          %max value for the amplitude
    
    elseif size(opt,2) == 5
        mY(k) = opt(4);          %min value for the time
        MY(k) = opt(5);          %max value for the time
        mXtmp=[];MXtmp=[];
        for i=1:Nsimu
           mXtmp=[mXtmp;s(i).mX(k)];
           MXtmp=[MXtmp;s(i).MX(k)];
        end
        mX(k)=min(mXtmp); MX(k)=mX(k)+dur;%max(MXtmp);
    
    else
        mYtmp=[];MYtmp=[];mXtmp=[];MXtmp=[];
        for i=1:Nsimu
           mYtmp=[mYtmp;s(i).mY(k)];
           MYtmp=[MYtmp;s(i).MY(k)];
           mXtmp=[mXtmp;s(i).mX(k)];
           MXtmp=[MXtmp;s(i).MX(k)];
        end
        mY(k)=min(mYtmp); MY(k)=max(MYtmp);
        mX(k)=min(mXtmp); MX(k)=mX(k)+dur;%max(MXtmp);
    end
end

%comment
disp(' ')
disp([num2str(nsS),' seismograms of Vn and Vp for ',num2str(Nsimu), ' type of simulations '])

%% PLOT SISMOGRAPHS 

for i=1:Nsimu

    iuc=0;
    
    %Sort points
    [A,sort1]=sort(s(i).x(s(i).ind));
    idx=s(i).ind(sort1);
    
    for cc=1:cf_uc:co
    iuc=iuc+1;

    %Sort points
    idxtmp=idx((iuc-1)*li+1:iuc*li);
    [B,sort2]=sort(s(i).z(idxtmp),'descend');
    idz=idxtmp(sort2);
    
    for j=1:li

        l=idz(j);
        
        %normalization elasticity
        if i==3
            s(i).vp(:,l)=s(i).vp(:,l)./max(abs(s(i).vp(:,l)))*max(abs(s(i-1).vp(:,l)));
            s(i).vn(:,l)=s(i).vn(:,l)./max(abs(s(i).vn(:,l)))*max(abs(s(i-1).vn(:,l)));
        end
        
        %find de dt difference between the a flat and a rough case
        if(i>2)
            idt2=find(abs(s(2).vn(:,l))==max(abs(s(2).vn(:,l))));
            idti=find(abs(s(i).vn(:,l))==max(abs(s(i).vn(:,l))));
            dtF=s(2).t(idt2)-s(i).t(idti);
        else
             dtF=0;
        end
        
        %plot
        subplot(li,co,co*(j-1)+cc:co*(j-1)+cc-1+cf_uc);     hold on
        plot(s(i).t+dtF,s(i).vn(:,l),'color',ctabvn(i,:),'Linewidth',1)
        plot(s(i).t+dtF,s(i).vp(:,l),'color',ctabvp(i,:),'Linewidth',1)

        %plot parameters
        k=find(l==s(i).ind);
        xlim([mX(k) MX(k)]);
        ylim([mY(k) MY(k)]);
        ax = gca; ax.FontSize = fs;
        if i==1, text(mX(k)+((MX(k)-mX(k))*0.02),MY(k)-((MY(k)-mY(k))*0.1),['z = ',num2str(round(10*s(i).z(l)/proZ)/10),' R_0'],'FontSize', fs);end
%         xlabel('time (s)','Fontsize',fs);
        if j< li; set(gca, 'xticklabel', []); else, xlabel('time (s)','Fontsize',fs); end
        if cc==1, ylabel('velocity (m/s)','Fontsize',fs); end
        if j == 1, if i==1, title(['x = ',num2str(round(10*uxc(iuc)/proZ)/10),' R_0'],'FontSize', fs),end, end
        if j == 1 && cc==1
            if i==1, text(mX(k)+((MX(k)-mX(k))*0.02),mY(k)+((MY(k)-mY(k))*(0.1*i)),['v_n',' damage'],'FontSize', fs,'color',ctabvn(i,:));end
            if i==1, text(mX(k)+((MX(k)-mX(k))*0.02),mY(k)+((MY(k)-mY(k))*(0.1*(i+1))),['v_p',' damage'],'FontSize', fs,'color',ctabvp(i,:));end
            if i==2, text(mX(k)+((MX(k)-mX(k))*0.2),mY(k)+((MY(k)-mY(k))*(0.1*(i-1))),['v_n',' elasticity'],'FontSize', fs,'color',ctabvn(i,:));end
            if i==2, text(mX(k)+((MX(k)-mX(k))*0.2),mY(k)+((MY(k)-mY(k))*(0.1*i)),['v_p',' elasticity'],'FontSize', fs,'color',ctabvp(i,:));end
        end

    end
    end

end

%save figures
if (save_tag>=1),export_fig(namefile,v_handle,'-eps');saveas(v_handle,[namefile,'.fig'],'fig');end
