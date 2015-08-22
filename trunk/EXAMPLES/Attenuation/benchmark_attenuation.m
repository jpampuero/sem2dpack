clear all
set(0,'DefaulttextFontSize',12)
   set(0,'DefaultaxesFontSize',12)
   set(0,'DefaultAxesLineWidth',2)
   set(0,'DefaultLineLineWidth',2)
   set(0,'DefaultLineMarkerSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbody=3;   

f0=18; fmin=1.8; fmax=180; 

QP=30;QS=20;

cs=2000;cp=3000;rho=2000;

%Source
F=-1;       %Amplitude of force       
t0=0.06;              % time shift

%Receiver
x=500;y=500;
r=sqrt(x^2+y^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nfrequency=2*Nbody-1;
w=zeros(Nfrequency,1);
wbody=zeros(Nbody,1);

AS=zeros(Nfrequency,Nbody);
QSInv=zeros(Nfrequency,1);
Y_beta=zeros(Nbody,1);

AP=zeros(Nfrequency,Nbody);
QPInv=zeros(Nfrequency,1);
Y_alpha=zeros(Nbody,1);

w0=2.0*pi*f0;
wmin=2.0*pi*fmin;
wmax=2.0*pi*fmax;

if (Nbody>1) 
    for i=1:Nfrequency
    w(i)=exp(log(wmin)+(i-1)*(log(wmax)-log(wmin))/(Nfrequency-1));
    end
else
    w(:)=w0;
end

for j=1:Nbody
     wbody(j)=w(2*j-1);
end 

for i=1:Nfrequency
    for j=1:Nbody
        AS(i,j)=(w(2*j-1)*w(i)+w(2*j-1)^2/QS)/(w(2*j-1)^2+w(i)^2);
        AP(i,j)=(w(2*j-1)*w(i)+w(2*j-1)^2/QP)/(w(2*j-1)^2+w(i)^2);
    end
end

QSInv(:)=1/QS;
QPInv(:)=1/QP;

Y_beta=pinv(AS)*QSInv;
Y_alpha=pinv(AP)*QPInv;

  RP1=1d0;
  RP2=0d0;
  RS1=1d0;
  RS2=0d0;

  for j=1:Nbody
     RP1=RP1-Y_alpha(j)/(1d0+(w0/wbody(j))^2);
     RP2=RP2+Y_alpha(j)*(w0/wbody(j))/(1d0+(w0/wbody(j))^2);
     RS1=RS1-Y_beta(j)/(1d0+(w0/wbody(j))^2);
     RS2=RS2+Y_beta(j)*(w0/wbody(j))/(1d0+(w0/wbody(j))^2);
  end 

  RP=sqrt(RP1^2+RP2^2);
  RS=sqrt(RS1^2+RS2^2);

  mu=rho*cs*cs;
  lambda=rho*(cp*cp-2d0*cs*cs);

  mu_inf=mu*(RS+RS1)/(2*RS^2);
  lambda_inf=(lambda+2d0*mu)*(RP+RP1)/(2*RP^2)-2d0*mu_inf;
  
  cs=sqrt(mu_inf/rho);
  cp=sqrt((lambda_inf+2*mu_inf)/rho);
  k_inf=lambda_inf+2/3*mu_inf;

wSeismax=2*pi*fmax;
wSeis=0:0.1:wSeismax;
wSeis=2*wSeis-wSeismax;  
mu_w=zeros(size(wSeis,2),1);
k_w=zeros(size(wSeis,2),1);
Y_k=(cp^2*Y_alpha-4/3*cs^2*Y_beta)/(cp^2-4/3*cs^2);
Y_mu=Y_beta;

for i=1:size(wSeis,2)
    mu_w(i)=mu_inf;
    k_w(i)=k_inf;
    for j=1:Nbody
        mu_w(i)=mu_w(i)-mu_inf*Y_mu(j)*wbody(j)/(wbody(j)+1i*wSeis(i));
        k_w(i)=k_w(i)-k_inf*Y_k(j)*wbody(j)/(wbody(j)+1i*wSeis(i));
    end
end

cp_w=sqrt((k_w+4/3*mu_w)/rho);
cs_w=sqrt(mu_w/rho);

u1=zeros(1,length(wSeis));
u2=zeros(1,length(wSeis));

for j=1:length(wSeis)
    if wSeis(j)>0
        G1=-1i*pi/2*(1/cp_w(j)^2*besselh(0,2,(wSeis(j)*r/(cp_w(j))))+ ...
            1./(wSeis(j)*r*cs_w(j))*besselh(1,2,(wSeis(j)*r/(cs_w(j))))- ...
            1./(wSeis(j)*r*cp_w(j))*besselh(1,2,(wSeis(j)*r/(cp_w(j)))));
        G2=1i*pi/2*(1/cs_w(j)^2*besselh(0,2,(wSeis(j)*r/(cs_w(j))))- ...
            1./(wSeis(j)*r*cs_w(j))*besselh(1,2,(wSeis(j)*r/(cs_w(j))))+ ...
            1./(wSeis(j)*r*cp_w(j))*besselh(1,2,(wSeis(j)*r/(cp_w(j)))));
        u1(j)=F/(2*pi*rho)*(x*y/r^2)*(G1+G2);
        u2(j)=F/(2*pi*rho)*(1/r^2)*(y^2*G1-x^2*G2);
    end
end

for j=1:length(wSeis)
        u1(j)=conj(u1(length(wSeis)+1-j));
        u2(j)=conj(u2(length(wSeis)+1-j));
end

S=sqrt(pi)*wSeis.^2/(4*(pi*f0)^3).*exp(-1i*wSeis*t0).*exp(-wSeis.^2/(4*(pi*f0)^2));

PHI_DX=u1.*S;
PHI_DY=u2.*S;

dt=1/(2*wSeismax)*2*pi;        

Sol_x=(ifft(fftshift(PHI_DX)))/dt;  
Sol_y=(ifft(fftshift(PHI_DY)))/dt;

t=(1:length(wSeis))*dt-dt ; 

for j=1:length(wSeis)/2
    t_res(j)=t(2*j-1);
end

for j=1:length(wSeis)/2
    SX_res(j)=Sol_x(2*j-1);
    SY_res(j)=Sol_y(2*j-1);
end

SX_res=abs(SX_res).*sign(real(SX_res));
SY_res=abs(SY_res).*sign(real(SY_res));

figure;
subplot(211);
plot(t_res,SX_res,'k');
title('(a) Horizontal displacement');
xlim([0.1 0.7]);
hold on
v=axis;
h=ylabel('Displacement (m)');
set(h,'Position',[v(1)-0.04 v(3)-3.5d-13 0]);
subplot(212);

plot(t_res,SY_res,'k');hold on
xlim([0.1 0.7]);

%%%%%%%% Results from SEM2DPACK %%%%%%%%

d=sem2d_read_seis;
subplot(211);
plot([1:d.nt]*d.dt,d.ux,'b');
xlim([0.1 0.7]);hold off
h1=legend('Analytical','SEM2DPACK');
set(h1,'fontsize',8);

subplot(212);
plot([1:d.nt]*d.dt,d.uz,'b');hold off
xlim([0.1 0.7]);
xlabel('Time (s)');
h1=legend('Analytical','SEM2DPACK');
set(h1,'fontsize',8);
title('(b) Vertical displacement');

