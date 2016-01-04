%plot figures showing biofilm deformation from a given position data set
%w49.91_f28_b0.4_cnd0.162_visc350_dt550_sigma0wSAR2_test12-16-15.mat
close
figure(1)
tstep=1735;
gplot3d(A,Xrstore(:,tstep*3+1:tstep*3+3));
hold on
scatter3(Xrstore(:,tstep*3+1),Xrstore(:,tstep*3+3),Xrstore(:,tstep*3+2),'k.');
hold off
axis 'equal'
axis([0 xlength -1.5*zlength 2.5*zlength 0 ylength]);
%h1=text(5,-1.5*zlength,ylength+2,'Top Plate','HorizontalAlignment','Left','Fontsize',12);
%set(h1,'rotation',15.5);
%h2=text(5,-1.5*zlength,2,'Bottom Plate','HorizontalAlignment','Left','Fontsize',12);
%set(h2,'rotation',15.5);


view(63,31);
xlabel('x ($\mu$m)','interpreter','latex','fontsize',14);
ylabel('z ($\mu$m)','interpreter','latex','fontsize',14);
zlabel('y ($\mu$m)','interpreter','latex','fontsize',14);

y=[ylength ylength ylength ylength]*0;
x=[0 xlength xlength 0];
z=[-1.5*zlength -1.5*zlength 2.5*zlength 2.5*zlength];
grid on
patch(x,z,y,'green','FaceAlpha',0.5);


y=[ylength ylength ylength ylength];
patch(x,z,y,'green','FaceAlpha',0.5);

title('Bacteria Positions at $t=0.063$','interpreter','latex','fontsize',14)
set(gca,'color','none');