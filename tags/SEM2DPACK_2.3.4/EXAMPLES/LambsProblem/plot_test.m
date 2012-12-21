% run first analyze_test.m

fac = 100; % factor to scale the difference

subplot(221)
plot(t,ux(:,1),t,uxa(:,1), t,(ux(:,1)-uxa(:,1))*fac,'--')
title('component X receiver 1')
legend('SEM2DPACK','EX2DDIR','diff \times 100','Location','SE')
subplot(222)
plot(t,ux(:,2),t,uxa(:,2), t,(ux(:,2)-uxa(:,2))*fac,'--')
title('component X receiver 2')
subplot(223)
plot(t,uz(:,1),t,uza(:,1),  t,(uz(:,1)-uza(:,1))*fac,'--')
title('component Z receiver 1')
subplot(224)
plot(t,uz(:,2),t,uza(:,2),  t,(uz(:,2)-uza(:,2))*fac,'--')
title('component Z receiver 2')
