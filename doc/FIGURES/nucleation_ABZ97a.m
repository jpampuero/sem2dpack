% plot describing the nucleation procedure kind = 2 in &BC_DYNFLT_TWF 
% expansion followed by contraction of a time-weakening region
% as in Andrews and Ben-Zion (JGR 1997, eqs 2)

subplot(3,1,[1 2])
L = 0.1;
x=[-0.6:0.005:0.6];
t=[0:0.005:1.2];
[T,X]=meshgrid(t,x);
r = t.*(1-t);
R = (abs(X)- repmat( r, length(x),1))/L;
MU = max( R, 0);
MU = min( MU, 1);
surf(X,T,MU)
view(2)
shading interp
axis([min(x) max(x) min(t) max(t)])
ylabel('Time / T')
title('Prescribed time-dependent friction coefficient')
%title('(\mu - \mu_d) /(\mu_s-\mu_d)')
% inverted colormap (to avoid spurious white edges)
map=colormap('gray');map=map(end:-1:1,:); colormap(map)
h=colorbar;
set(h,'YTick',[0 1],'YTickLabel',{'mu_d' 'mu_s'});
% sorry, no easy way to use LaTeX in labels ! 
% (but see format_labels and tick2tex in Matlab Central File Exchange)

k=floor(length(t)/2);
line([min(x) max(x)], [1 1 ]*t(k), 'LineStyle','--','Color','b');

subplot(3,1,3)
plot(x,MU(:,k))
%axis([min(x) max(x) -0.2 1.2])
axis([min(x) 0.8 -0.2 1.2])
xlabel('Position along fault / (V x T)')
set(gca,'YTick',[0 1],'YTickLabel',{'mu_d' 'mu_s'});
