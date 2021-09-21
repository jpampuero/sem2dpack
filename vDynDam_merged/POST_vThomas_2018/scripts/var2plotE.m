% Script to download data 

% Marion Thomas last modified september 2018

%CALLS: 

%==========================================================================
%% 
%Download damage variable
elm = sem2d_snapshot_read('elm',oo(o(i)),datadir);
elm0= sem2d_snapshot_read('elm',oo(2),datadir);

%default values for plot options
barvar='';
vname='';
writcol=wcol;
colorL=[0 0 0];

%Find the variable to plot
if strcmp('ax',var) == 1 || strcmp('az',var) == 1
    v = sem2d_snapshot_read(var(i,:),oo(o),datadir);
    barvar='m^2/s';vname=var;
    
elseif strcmp('vx',var) == 1
    v = sem2d_snapshot_read(var(i,:),oo(o),datadir);
    barvar='m/s';vname='v_p';

elseif strcmp('vz',var) == 1
    v = sem2d_snapshot_read(var(i,:),oo(o),datadir);
    barvar='m/s';vname='v_n';
    
elseif strcmp('nv',var) == 1
    v1 = sem2d_snapshot_read('vz',oo(o(i)),datadir);
    v2 = sem2d_snapshot_read('vx',oo(o(i)),datadir);
    v=(v1.^2+v2.^2).^(1/2);
    barvar='m/s';vname='\sqrt{v_n^2+v_p^2}';
    
elseif strcmp('s11',var) == 1 || strcmp('s12',var) == 1 || strcmp('s22',var) == 1
    v = sem2d_snapshot_read(var(i,:),oo(o),datadir);
	v=v/1e6; barvar='MPa';vname=var;
    
elseif strcmp('e11',var) == 1 || strcmp('e12',var) == 1 || strcmp('e22',var) == 1
    v = sem2d_snapshot_read(var(i,:),oo(o),datadir);
    barvar='';vname=var;

elseif strcmp('invI',var) == 1
    v10=elm0.invI;mI0=max(max(max(v10)));
    v=(elm.invI-mI0)/abs(mI0);
    barvar='';vname='\Delta I_1 = (I_1-I_1^{0})/|I_1^{0}|';
    
elseif strcmp('invII',var) == 1
    v20=elm0.invII;mII0=abs(max(max(max(v20))));
    v=(elm.invII-mII0)/abs(mII0);
    barvar='';vname='\Delta J_2 = (J_2-J_2^{0})/|J_2^{0}|';
  
elseif strcmp('CFF',var) == 1
    v=(elm.invII+elm.invI*0.6)/1e6;
    barvar='MPa';vname='\Delta CF = J_2+0.6*I_1';

else
    disp([var ' is not a valid name for an output variable'])
end