ngll=5;
nex=20;
zmin=-1e3;
xmin=-1e3; xmax=1e3;

ns=nex*(ngll-1)+1;
clear table
table(:,1) = xmin+(xmax-xmin)/(ns-1)*(0:ns-1)';
table(:,2) = zmin;
table(:,3) = (table(:,1)-xmin)/4e3;  % = x / horizontal_speed

save('input_sources.tab', 'table', '-ascii')
