FMAX = 20;	% maximum frequency to be resolved
NEPW = 1; 	% minimum number of elements per wavelength
XLIM = [-50 50];	% horizontal limits of the model
NBOT = 1; 	% bottom layer thickness = NBOT * minimum wavelength
                % free surface is z=0

[model.rho, model.cp, model.cs, model.h]= ...
  textread('RomaNorte_model.tab','%*s%f%f%f%*f%f', 'commentstyle','matlab');

model.h(end) = NBOT*model.cs(end)/FMAX;
model.xlim = XLIM;

p = set_mesh_layers(model,FMAX,NEPW);

