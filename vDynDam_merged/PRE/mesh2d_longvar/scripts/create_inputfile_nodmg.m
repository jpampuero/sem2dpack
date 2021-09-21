% Script to write the input file Par.inp for SEM2DPACK, for a pure elastic
% case.

% Marion Thomas, last modified June 2018

%CALLS: 

%==========================================================================
%% Create file

namefile='Par.inp';
wfile = [rootD,namefoldE,'/',namefile];
fid = fopen(wfile,'wt');

%% Headline

HEAD{1}=['# This is an input file for SEM2DPACK, using the mesh: ', namemesh,'.mesh2d'];
fprintf(fid,[HEAD{1},' \n']);%; fprintf(fid,' \n');
fprintf(fid,' \n');

%% #----- General parameters ----------------

fprintf(fid,'#----- General parameters ---------------- \n');

fprintf(fid,'&GENERAL iexec=1'); %iexec=0 corresponds to the check mode
fprintf(fid,[', ngll=',num2str(ngll)]);
fprintf(fid,[', fmax=',num2str(fmax)]);
fprintf(fid,[', ndof=',num2str(ndof),' ,']);
fprintf(fid,' \n');

fprintf(fid,'  title = ''Slip-weakening dynamic rupture with off-fault damage''');
fprintf(fid,', verbose=''1011''');
fprintf(fid,', ItInfo = 1 / \n');
fprintf(fid,' \n');

%% #----- Build the mesh ---------------------------

fprintf(fid,'#----- Build the mesh --------------------------- \n');


fprintf(fid,'&MESH_DEF  method = ''MESH2D'' / \n');
fprintf(fid,['&MESH_MESH2D file=''',namemesh,'.mesh2d'' / \n']);
fprintf(fid,' \n');

%% #---- Material parameters --------------

%Find the number of different tag
Ndom=length(domain);
tagval=[];
for i=1:Ndom
    tagval=[tagval;unique(domain(i).etag)];
end
Ntag=numel(unique(tagval));


fprintf(fid,'#---- Material parameters -------------- \n');
j=1;
for i=1:Ntag
    
if Nmat==2 && D0z(i) >= 0, j=2; end

fprintf(fid,['&MATERIAL tag=',num2str(i),', kind=''ELASM'' / \n']);

fprintf(fid,['&MAT_ELASTICM rho=',num2str(rho(j))]);
fprintf(fid,[', cs=',num2str(cs(j))]);
fprintf(fid,[', cp=',num2str(cp(j)),' \n']);
fprintf(fid,'           ');
fprintf(fid,['e0=',num2str(e0(1,j)),', ',num2str(e0(2,j)),', ',num2str(e0(3,j)),' / \n']);
fprintf(fid,' \n');

end

%% #----- Boundary conditions ---------------------

fprintf(fid,['#----- Boundary conditions --------------------- \n']);

fprintf(fid,'&BC_DEF  tags = 5,6 , kind = ''DYNFLT'' / \n');
fprintf(fid,['&BC_DYNFLT friction=''SWF'',Szz=',num2str(round(Szz)),...
    ',Sxz=',num2str(round(Sxz)),',Sxx=',num2str(round(Sxx)),',Tn=0.0d7,TtH=''PWCONR'' / \n']);
fprintf(fid,['&DIST_PWCONR num=2, ref=',num2str(round(nuclocX)),',',num2str(round(nuclocZ)), '/ # Initial shear stress \n']);
fprintf(fid,['   ',num2str(nucsize),'\n']);
fprintf(fid,[num2str(round(dTt_nuc)),'  ','0.0d7','\n']);

fprintf(fid,['&BC_DYNFLT_SWF Dc=',num2str(Dc)]);
fprintf(fid,[', MuS=',num2str(MuS)]);
fprintf(fid,[', MuD=',num2str(MuD),' / \n']);

fprintf(fid,['&BC_DYNFLT_NOR kind=',num2str(kind)]);
fprintf(fid,[', T=',num2str(T),' / \n']);
fprintf(fid,' \n');

fprintf(fid,'&BC_DEF  tag = 1 , kind = ''ABSORB'' / \n');
fprintf(fid,'&BC_DEF  tag = 2 , kind = ''ABSORB'' / \n');
fprintf(fid,'&BC_DEF  tag = 3 , kind = ''ABSORB'' / \n');
fprintf(fid,'&BC_DEF  tag = 4 , kind = ''ABSORB'' / \n');

fprintf(fid,' \n');
 
%% #---- Time scheme settings ----------------------

fprintf(fid,['#---- Time scheme settings ---------------------- \n']);

% &TIME  TotalTime=5.0d0, courant = 0.55d0, kind='leapfrog' /  ! T=30
fprintf(fid,['&TIME  TotalTime=',num2str(TOTtime)]);
fprintf(fid,', courant = 0.55d0, kind=''leapfrog'' / \n');
fprintf(fid,' \n');
 
%% #----- Receivers ---------------------------------

fprintf(fid,['#----- Receivers --------------------------------- \n']);

fprintf(fid,['&REC_LINE  file=''',sensorfile,'''']);
fprintf(fid,[', isamp=',num2str(isamp)]);
fprintf(fid,[', field=''',fieldSen,''' / \n']);
fprintf(fid,' \n');

%% #--------- Plots settings ----------------------
% &SNAP_DEF itd=400, fields ='ESVA', ps=F /
fprintf(fid,['#--------- Plots settings ---------------------- \n']);

fprintf(fid,['&SNAP_DEF itd=',num2str(itd)]);
fprintf(fid,[', fields=''',fieldSnap,'''']);
fprintf(fid,[', ps=',ps,' / \n']);
fprintf(fid,' \n');

fclose(fid);


