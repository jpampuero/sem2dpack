% Load the finite element data
Knods = load('ElmtNodes_sem2d.tab');
Coorg = load('MeshNodesCoord_sem2d.tab');
NElem = size(Knods,1);
Ex=reshape(Coorg(Knods,1),NElem,4)'; % +1.162e6; 
Ey=reshape(Coorg(Knods,2),NElem,4)';

% Load the check data 
Reso  = load('Resolution_sem2d.tab');
Stab  = load('Stability_sem2d.tab');

% Plots
subplot(2,3,[1 2]);fill(Ex,Ey,1./Stab'); axis equal; colorbar
title('Stability: 1/CFL \propto local \Deltat')
subplot(2,3,3), hist(1./Stab)
subplot(2,3,[4 5]);fill(Ex,Ey,Reso'); axis equal; colorbar
title('Resolution: number of nodes per \lambda_{min}')
subplot(2,3,6), hist(Reso)
