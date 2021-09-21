% Script to download data 

% Marion Thomas last modified september 2018

%CALLS: 

%==========================================================================
%% 
%Download damage variable
dyD = sem2d_snapshot_read('dyD',oo(o(i)),datadir);
dyD0= sem2d_snapshot_read('dyD',oo(2),datadir);

%default values for plot options
barvar='';
vname='';
writcol=wcol;
colorL=[0 0 0];

%Find the variable to plot
if strcmp('D',var)==1
    v=dyD.D; 
    %v(v<=0.25)=0;NaN;
    barvar='D';vname='Damage';
    colorL=[1 0 0];
    
elseif strcmp('cs',var) == 1
    v=dyD.cs;
    cs0=max(max(max(v)));
    v=100*(cs0-v)/cs0;
    barvar='%';vname='\%\ of\ reduction\ in\ c_s';

elseif strcmp('cp',var) == 1
    v=dyD.cp;
    cp0=max(max(max(v)));
    v=100*(cp0-v)/cp0; 
    barvar='%';vname='\%\ of\ reduction\ in\ c_p';
    
elseif strcmp('dldt',var) == 1
    v=dyD.dldt;
    
elseif strcmp('R',var) == 1
    maxV=3.5;minV=0.5;
    v=dyD.R;
    barvar='R';vname='Regime';
%     writcol='white';colorL=[1 0 0];%colorL=[1 1 1];
    
elseif strcmp('ax',var) == 1 || strcmp('az',var) == 1
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
    v10=dyD0.invI;mI0=max(max(max(v10)));
    v=(dyD.invI-mI0)/abs(mI0);
%     barvar='';vname='(\sigma-\sigma_0)/|\sigma_0|';
    barvar='';vname='\Delta I_1 = (I_1-I_1^{0})/|I_1^{0}|';
    
elseif strcmp('invII',var) == 1
    v20=dyD0.invII;mII0=abs(max(max(max(v20))));
    v=(dyD.invII-mII0)/abs(mII0);
%     barvar='';vname='(\tau-\tau_0)/\tau_0';
    barvar='';vname='\Delta J_2 = (J_2-J_2^{0})/|J_2^{0}|';

elseif strcmp('CFF',var) == 1
    v=(dyD.invII+dyD.invI*0.6)/1e6;
%     barvar='MPa';vname='(\tau+0.6*\sigma)';
    barvar='MPa';vname='\Delta CF = J_2+0.6*I_1';

elseif strcmp('er11',var) == 1
    v=dyD.er11;
    barvar='.s-1';vname=var;
    
elseif strcmp('er22',var) == 1
    v=dyD.er22;
    barvar='.s-1';vname=var;
    
elseif strcmp('er12',var) == 1
    v=dyD.er12;
    barvar='.s-1';vname=var;
    
else
    disp([var ' is not a valid name for an output variable'])
end