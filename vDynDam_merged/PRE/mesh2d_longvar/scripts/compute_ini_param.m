% Script to compute initial strain based on fault angle with principal
% stresses as well as the others necessary input parameters

% Marion Thomas, last modified August 2018

%CALLS: 

%==========================================================================
%% DYNAMIC DAMAGE PARAMETERS


%Own definition of D0
if own_tag==1;
    D0p=ownD0;
    D0pz=ownD0p; 
    %correction of D0pz based on the resolution
    D0pz=resocorrec(D0pz,h,Q);
end

%initial Damage
if numel(D0p)>0
    D0ztemp=[BOX(3);D0pz;BOX(4)];
    D0temp=[D0bot;D0p;D0top];
    [z,iz]=sort(D0ztemp);
    D0=D0temp(iz);
    D0z=D0ztemp(iz);
elseif D0bot~=D0top
    D0=[D0bot;D0top];
    D0z=[BOX(3);BOX(4)];
elseif numel(rho)==2 && D0bot==D0top
    D0=[D0bot;D0top];
    D0z=[BOX(3);BOX(4)];
else
    D0=D0bot;
    D0z=[BOX(3);BOX(4)];
end

% %Correction based on the resolution
% D0z=resocorrec(D0z,h,Q);

%others parameters
lambda=rho.*(cp.*cp - 2d0.*cs.*cs);
mu=rho.*cs.*cs;
Phi = 0.5*atand(1/fs);% angle that the penny shape cracks make to the largest most compressive stress sigmal
alpha = cosd(Phi);  % projection of the crack radius in a vertical plane parallel to the direction of sigma1
Nv=D0*(3/4)/pi*(alpha*a0)^(-3); % Number of cracks
nu=lambda./(2*(lambda+mu));
cr = ((0.87 + 1.12.*nu)./(1+nu)).*cs;

%number of material
Nmat=ceil((numel(rho)+numel(cs)+numel(cp))/3);
 
%% DETERMINE INITIAL STRAIN, SXZ, SXX AND MU0

kronecker = eye(size(zeros(3,3)));

%Mohr circle but for compression negative convention
r=(PS(3)-PS(1))/2;
O=(PS(3)+PS(1))/2;

%Background stress
Bck_stress(1)= O + r*cosd(2*ang);
Bck_stress(2)= O - r*cosd(2*ang);
Bck_stress(3)= r*sind(2*ang);


%PLANE STRAIN (like in the code)

%declare variables
e0=zeros(3,Nmat);
e=zeros(3,3,Nmat);
S=zeros(3,3,Nmat);
tmpe=zeros(3,3);

for i=1:Nmat
    
    %In term of strain
    sigma = Bck_stress(1)+Bck_stress(2);
    e0(1,i)= Bck_stress(1)/(2*mu(i))-(nu(i)/(2*mu(i)*(1+nu(i))))*sigma;
    e0(2,i)= Bck_stress(2)/(2*mu(i))-(nu(i)/(2*mu(i)*(1+nu(i))))*sigma;
    e0(3,i)= Bck_stress(3)/(2*mu(i));

    %Strain field based on e0 (no e22)
    e(1,1,i)=e0(1);
    e(3,3,i)=e0(2);
    e(1,3,i)=e0(3);
    e(3,1,i)=e0(3);

    %Compute stress
    tmpe=e(:,:,i);
    [~,ds] = eigs(tmpe);
    eps = e(1,1,i)+e(2,2,i)+e(3,3,i);trace(ds);
    S(:,:,i)=2*mu(i)*tmpe+lambda(i)*eps*kronecker;
    Mu0=abs(S(1,3,i)/S(3,3,i));

end

%% PAREMETERS TO MODIFY IN PAR.INP

disp('Parameters for the Par.inp file')
disp(' ')

disp('#---- Material parameters --------------')

% fs=fs*1 
disp(['nu = ',num2str(nu)])
disp(['cr = ',num2str(cr)])
Nv=Nv*1;
disp(['Nv = ',num2str(Nv')])
disp(['e0 = ',num2str(e0(:)')])
disp(' ')

disp('#---- Boundary conditions --------------')
Szz=S(3,3); disp(['Szz = ',num2str(Szz./1e6)])
Sxz=S(1,3); disp(['Sxz = ',num2str(Sxz./1e6)])
Sxx=S(1,1); disp(['Sxx = ',num2str(Sxx./1e6)])
Mu0=Mu0*1;  disp(['Mu0 = ',num2str(Mu0)])
if (Mu0<MuD)
    disp('pb Mu0 has to be greater than MuD')
end
if (Mu0>MuS)
    disp('pb MuS has to be greater than Mu0')
end 
disp(' ')

disp('#---- fault properties -----------------')
tauS=MuS*-Szz;      disp(['tauS = ',num2str(tauS./1e6)])
tauD=MuD*-Szz;      disp(['tauD = ',num2str(tauD./1e6)])
tau0=Sxz;           disp(['tau0 = ',num2str(tau0./1e6)])
Strength=(tauS-tau0)/(tau0-tauD);       disp(['S parameter = ',num2str(Strength)])
disp(' ')

%% Resolution

dh=h/(ngll^1);

if (ndof==2)
    mu_s=mu./(1-nu);
else
    mu_s=mu;
end

pro_zone = (9*pi/32)*mu_s/((-MuS+MuD)*S(3,3));

disp('#---- Resolution -----------------------')
disp(['process zone (R0) = ', num2str(pro_zone)])
disp(['h = ', num2str(h)])
disp(['grid size = ',num2str(dh)])
disp(' ')


%% stress invariant

disp('#---- Stress invariant -----------------')
% S(2,2)=0;

for k=1:Nmat
    sig=trace(S(:,:,k))./3;

    delta=eye(3,3);
    devS=S-(sig*delta);
    dev=zeros(3,3);
    for i=1:3
        for j=1:3
            dev(i,j)=devS(i,j)*devS(j,i);
        end
    end
    tau = sqrt(0.5*sum(sum(dev))) ;

    disp(['sigma = ', num2str(sig/1e6)])
    disp(['tau = ', num2str(tau/1e6)])
    disp(' ')
end

