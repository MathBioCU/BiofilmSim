load('../../ComplianceData.mat');

figure(1)
P1=loglog(Time,J0_5,'-b');
hold on
P2=loglog(Time,J0_1,'-g');

h=legend([P1,P2],'$J(t;0.5)$ Experimental','$J(t;0.1)$ Experimental');
set(h,'interpreter','latex','location','southeast');