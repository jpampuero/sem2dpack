% ************************************************************************
% This function plot the fourrier analysis of the particules velocity ( cf
% plot_velocity) for the different "seismometers" (point outputs, definied in the input file)

% Marion Thomas, last modified October 2019

% CALLS: FAS.m; export_fig.m;

%% ************************************************************************
function [] = plot_FAS_ratio_VpVn(Pcoord,data,fault,tend,save_tag,namef,opt)

%% INPUTS VARIABLES

%Pcoord:	seismographs locations
%data:      seismograms (potentially including several simulations)
%save_tag:  saving option (1=save, 0=just plot, 2=save just .fig)
%namef      field to be plot
%field      definie which component (z or x) should be plotted

%options
proZ = opt(2);	%process zone size
fs = opt(1);	%font size
fcutoff = opt(3);      %Frequency cut off
taperlengthpercent = opt(4); %Cosine Taper

%% PLOT PARAMETERS I

%Launch figure
v_handle =figure('Position',[0 500 1600 1000],'PaperOrientation','landscape');  
fv=v_handle;
set(fv,'PaperOrientation','landscape');
set(fv,'PaperUnits','normalized');
set(fv,'PaperPosition', [0 0 1 1]);
namefile=[namef,'FASratio_vpvn'];

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
ctabvn=[0 0 0;1 0 0;0.4 0 0;0 0 0];
% ctabvp=[0 0 1;0 0 0.4;0 0 0];
ctabvp=[0 0 0;0.3 0.6 1;0 0 1;0 0 0];

%% DOWNLOAD DATA

%number of simulations
Nsimu=size(data,1);
s(1).nt=0;
s(1).ns=0;

%Number of ratio
Nratio=factorial(Nsimu)/(factorial(2)*factorial(Nsimu-2));

%Ration combination
ratioM=zeros(Nratio,2);
r=0;
for i=1:Nratio
    for j=i+1:Nsimu
        r=r+1;
        ratioM(r,:)=[i,j];
    end
end

for i=1:Nsimu
    
%field to compute the FAS
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
    end
    
    %for the FAS    
	iMX=find(abs(s(i).t-s(i).tend)== min(abs(s(i).t-s(i).tend)));
    s(i).srate = 1/max(diff(s(i).t));   %sampling rate 
    s(i).taperstart =s(i).t(iMX)*(100-taperlengthpercent)/100;
    s(i).taperdur = s(i).t(iMX)-s(i).taperstart;


end

%% PLOT PARAMETERS II

%min/max for amplitude and time
mY = opt(5);          %min value for the time
MY = opt(6);          %max value for the time
mX = -2;
MX = log10(fcutoff);

%comment
disp(' ')
disp([num2str(nsS),' FAS analaysis of Vn and Vp for seismograms for ',num2str(Nsimu), ' type of simulations '])


%% Compute the FAS

for i=1:Nsimu

    iuc=0;
    
    %Sort points
    [A,sort1]=sort(s(i).x(s(i).ind));
    idx=s(i).ind(sort1);
    
    %for FAS
	[a,b] = butter(1,fcutoff*0.5/s(i).srate);

    for cc=1:cf_uc:co
    iuc=iuc+1;

    %Sort points
    idxtmp=idx((iuc-1)*li+1:iuc*li);
    [B,sort2]=sort(s(i).z(idxtmp),'descend');
    idz=idxtmp(sort2);
    
    for j=1:li

        l=idz(j);
        
        %FAS
        tempvn = filter(a,b,s(i).vn(:,l));
        tempvp = filter(a,b,s(i).vp(:,l));
        [s(i).fasvn(:,l),s(i).f1vn(:,l),v2] = FAS(s(i).t,tempvn, s(i).taperstart, s(i).taperdur);    
        [s(i).fasvp(:,l),s(i).f1vp(:,l),v2] = FAS(s(i).t,tempvp, s(i).taperstart, s(i).taperdur);    
    end
    end

end

%% PLOT the ratio of FAS

numel_dataf1=[];
for ns=1:Nsimu
    numel_dataf1=[numel_dataf1;numel(s(ns).f1vn(:,l))];
end
mini=find(numel_dataf1 == min(numel_dataf1),1);


for k=2:Nratio

    iuc=0;
    i1=ratioM(k,1);
    i2=ratioM(k,2);
    
    %Sort points
    [A,sort1]=sort(s(i).x(s(i1).ind));
    idx=s(i1).ind(sort1);

	co2=co*2;
    
    for cc=1:cf_uc:co
    iuc=iuc+1;

    %Sort points
    idxtmp=idx((iuc-1)*li+1:iuc*li);
    [B,sort2]=sort(s(i).z(idxtmp),'descend');
    idz=idxtmp(sort2);
    
    for j=1:li

        l=idz(j);
        
        %resample and normalized for Vn
        xsamp=s(ratioM(k,mini)).f1vn(:,l);
        ysamp1=interp1(s(i1).f1vn(:,l),s(i1).fasvn(:,l)./s(i1).fasvn(1,l),xsamp);
        ysamp2=interp1(s(i2).f1vn(:,l),s(i2).fasvn(:,l)./s(i2).fasvn(1,l),xsamp);

        %plot
        subplot(li,co2,co2*(j-1)+cc:co2*(j-1)+cc-1+cf_uc);     hold on
        plot(log10(xsamp),log10(ysamp1./ysamp2),'-','color',ctabvn(k,:),'Linewidth',1)
       
        %plot parameters
        xlim([mX MX]); ylim([mY MY])
        ax = gca; ax.FontSize = fs;
        if k==3; text(mX+((MX-mX)*0.02),MY-((MY-mY)*0.1),['z = ',num2str(round(10*s(i).z(l)/proZ)/10),' R_0'],'FontSize', fs);end
        text(mX+((MX-mX)*0.02),MY-((MY-mY)*0.1*(k+1)),['ratio d',num2str(i1),'/d',num2str(i2)],'FontSize', fs,'color',ctabvn(k,:));
        if j< li; set(gca, 'xticklabel', []); else, xlabel('log f','Fontsize',fs); end
        if cc>cf_uc, set(gca, 'yticklabel', []); else, ylabel(['log of velocity' ],'Fontsize',fs),end
        if j == 1, title(['x = ',num2str(round(10*uxc(iuc)/proZ)/10),' R_0'],'FontSize', fs),end

    
        %resample and normalized for Vp
        mini=find([numel(s(i1).f1vp(:,l)),numel(s(i2).f1vp(:,l))] == min([numel(s(i1).f1vp(:,l)),numel(s(i2).f1vp(:,l))]),1);
        xsamp=s(ratioM(k,mini)).f1vp(:,l);
        ysamp1=interp1(s(i1).f1vp(:,l),s(i1).fasvp(:,l)./s(i1).fasvp(1,l),xsamp);
        ysamp2=interp1(s(i2).f1vp(:,l),s(i2).fasvp(:,l)./s(i2).fasvp(1,l),xsamp);

        %plot
        subplot(li,co2,co2*(j-1)+cc+co:co2*(j-1)+cc-1+cf_uc+co);     hold on
        plot(log10(xsamp),log10(ysamp1./ysamp2),'-','color',ctabvp(k,:),'Linewidth',1)
       
        %plot parameters
        xlim([mX MX]); ylim([mY MY])
        ax = gca; ax.FontSize = fs;
        if k==3; text(mX+((MX-mX)*0.02),MY-((MY-mY)*0.1),['z = ',num2str(round(10*s(i).z(l)/proZ)/10),' R_0'],'FontSize', fs);end
        text(mX+((MX-mX)*0.02),MY-((MY-mY)*0.1*(k+1)),['ratio d',num2str(i1),'/d',num2str(i2)],'FontSize', fs,'color',ctabvp(k,:));
        if j< li; set(gca, 'xticklabel', []); else, xlabel('log f','Fontsize',fs); end
        set(gca, 'yticklabel', []);
        if j == 1, title(['x = ',num2str(round(10*uxc(iuc)/proZ)/10),' R_0'],'FontSize', fs),end

    
    end
    end

end

%save figures
if (save_tag>=1),export_fig(namefile,v_handle,'-eps');saveas(v_handle,[namefile,'.fig'],'fig');end