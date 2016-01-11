load('../../Miscellaneous/ComplianceData.mat');

figure(1)
P2=loglog(Time,J0_2,'-.b');
hold on
P1=loglog(Time,J0_1,'-.g');

h=legend([P1,P2],'$J(t;0.1)$ Experimental','$J(t;0.2)$ Experimental');
set(h,'interpreter','latex','location','southeast');