model_dir = ['../EXAMPLES/2.5D_inplane/'];

g = sem2d_read_specgrid(model_dir);

sem2d_snapshot_movie('vz',g,'DATADIR',model_dir,'TPAUSE',1 )
