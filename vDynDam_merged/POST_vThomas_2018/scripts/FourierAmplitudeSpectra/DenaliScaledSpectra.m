clear all

load DenaliLab.mat

% Cosine Taper
taperlengthpercent = 10;

% Resample lab data at the same sampling frequency as the Denali data
% Has been checked and this level of coarsening does not affect results too
% much

tend = 39;
index = find(denali.t>=tend,1);

tsamp = denali.t(1:index);

vpden = interp1(denali.t,denali.vp,tsamp);
vnden = interp1(denali.t,denali.vn,tsamp);
vzden = interp1(denali.t,denali.vz,tsamp);


vpsamp = interp1(scaledSS.tp,scaledSS.vp,tsamp);
vnsamp = interp1(scaledSS.tn,scaledSS.vn,tsamp);
vzsamp = interp1(scaledSS.tz,scaledSS.vz,tsamp);

samplingrate = 1/max(diff(tsamp));

taperstart = tsamp(end)*(100-taperlengthpercent)/100;
taperdur = tsamp(end)-taperstart;

% var = {'sqrt(vzsamp.^2+vnsamp.^2+vpsamp.^2)' 'sqrt(denali.vp.^2+denali.vn.^2+denali.vz.^2)'};
% var = {'vnsamp','denali.vn'};
var = {'vzsamp','denali.vz'};
% var = {'vpsamp','denali.vp'};

clr = {'b','r','g'};
clf;
fcutoff = 20;
[a,b] = butter(1,fcutoff*0.5/samplingrate);


for i=1:2
    eval(['temp = ' var{i} ';']);
    
%     out = fourier(taper(temp,taperlengthpercent),samplingrate);
%     f1 = out(:,1);
%     fas1 = out(:,2);

    temp = filter(a,b,temp);
    
    [fas1,f1,v2] = FAS(tsamp,temp, taperstart, taperdur);    

    hold all
    
    plot(log10(f1),log10(fas1),clr{i})
    xlabel('log$_{10} \ f$');
    ylabel('log$_{10} \ $FAS');

%     i1 = find(log10(f1)>log10(3),1);
%     i2 = find(log10(f1)>log10(30),1);
%     x = log10(f1);
%     y = log10(fas1);
%     
%     xband = x(i1:i2);
%     yband = y(i1:i2);
%     P = polyfit(xband,yband,1);
%     plot(xband,P(1)*xband+P(2),'--k','Linewidth',2);
%     
%     slopex = xband(round(0.5*(i2-i1)));
%     slopey = P(1)*slopex+P(2);
%     
%     if i==1
%         text(slopex,slopey-1,['$f^{' num2str(P(1),'%3.2f') '}$']);
%     else
%         text(slopex,slopey+1,['$f^{' num2str(P(1),'%3.2f') '}$']);
%     end
% 
%     i1 = find(log10(f1)>log10(0.01),1);
%     i2 = find(log10(f1)>log10(3),1);
%     x = log10(f1);
%     y = log10(fas1);
%     
%     xband = x(i1:i2);
%     yband = y(i1:i2);
%     P = polyfit(xband,yband,1);
%     plot(xband,P(1)*xband+P(2),'--g','Linewidth',2);
% 
%     slopex = xband(round(0.5*(i2-i1)));
%     slopey = P(1)*slopex+P(2);
%     
%     if i==1
%         text(slopex,slopey-1,['$f^{' num2str(P(1),'%3.2f') '}$']);
%     else
%         text(slopex,slopey+1,['$f^{' num2str(P(1),'%3.2f') '}$']);
%     end
    xlim([-2 log10(fcutoff)]);
    ylim([-3 3]);
end

legend('Lab','Denali')

%%
% clf;
% tsamp = denali.t;
% vsamp = denali.vn;
% vsamp = taper(vsamp,taperlengthpercent);
% samplingrate = 1/max(diff(tsamp));
% 
% window = 500;
% noverlap = 40;
% nfft = max(256,nextpow2(samplingrate/window));
% [S,F,T,P] = spectrogram(vsamp,window,noverlap,nfft,samplingrate);
% subplot(2,1,1)
% hold all
% plot(tsamp,vsamp*50+15,'-k')
% surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
% view(0,90);
% xlabel('Time (Seconds)'); ylabel('Hz');
% 
% tsamp = denali.t;
% vsamp = interp1(scaledSS.tp,scaledSS.vn,tsamp);
% vsamp = taper(vsamp,taperlengthpercent);
% samplingrate = 1/max(diff(tsamp));
% 
% nfft = max(256,nextpow2(samplingrate/window));
% [S,F,T,P] = spectrogram(vsamp,window,noverlap,nfft,samplingrate);
% subplot(2,1,2)
% hold all
% plot(tsamp,vsamp*50+15,'-k')
% surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
% view(0,90);
% xlabel('Time (Seconds)'); ylabel('Hz');


