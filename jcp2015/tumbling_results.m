%tumbling rotation pca for ellipsoid

t=-1:h:1;

for i=1:5:c3
Z=Xrstore(:,3*i-2:3*i)-ones(length(X),1)*mean(Xrstore(:,3*i-2:3*i));
Zr=Z*pca(Z);
V=pca(Z);

Xs=Xrstore(:,3*i-2:3*i);
%Xs=Zr;
mX=mean(Xs);
scatter3(Xs(:,1),Xs(:,3),Xs(:,2),'k.')
hold on; V=V';
plot3(V(1,1)*t+mX(1),V(1,3)*t+mX(3),V(1,2)*t+mX(2),'r','linewidth',2)
plot3(V(2,1)*t+mX(1),V(2,3)*t+mX(3),V(2,2)*t+mX(2),'g','linewidth',2)
plot3(V(3,1)*t+mX(1),V(3,3)*t+mX(3),V(3,2)*t+mX(2),'b','linewidth',2)

plot3(t+mX(1),0*t+mX(3),0*t+mX(2),'g','linewidth',1);
plot3(0*t+mX(1),t+mX(3),0*t+mX(2),'g','linewidth',1);
plot3(0*t+mX(1),0*t+mX(3),t+mX(2),'g','linewidth',1);

[x,y,z]=ellipsoid(0,0,0,...
    2*std(Xs(:,1)),2*std(Xs(:,3)),2*std(Xs(:,2)),20);

% rotate to eigenframe
for j=1:length(z)
    for k=1:length(z)
        G=[x(j,k); y(j,k); z(j,k)];
        G=V'*G;
        x(j,k)=G(1);
        y(j,k)=G(2);
        z(j,k)=G(3);
    end
end

x=x+mX(1);
y=y+mX(2);
z=z+mX(3);

AA=surf(x,z,y,'facecolor','g','edgecolor','none');
alpha(AA,0.5);
axis('equal')
 %axis([-xlength 3*xlength -100*zlength 1*zlength -ylength ylength]); %Xs=X()
%axis([-xlength/2 xlength/2 -zlength/2 zlength/2 -ylength/4 ylength/4]); %Xs=Zr
camlight

view(-14,38);
view(-90,10);
str1=num2str(round(std(Xs(:,1))*1000)/1000);
str2=num2str(round(std(Xs(:,2))*1000)/1000);
str3=num2str(round(std(Xs(:,3))*1000)/1000);
str=['$\sigma_x=$', str1, '  $\sigma_y=$ ',str2, '  $\sigma_z=$ ', str3];

title(str,'Interpreter','Latex','Fontsize',16)
shg
pause(0.1)
hold off
end