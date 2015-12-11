subplot(2,1,1)
plot(fStrain(1:c3))
str = ['Strain, \omega=', num2str(w)];
title(str);
hold on
[mstrain,tstrain]=max(fStrain);
plot(tstrain*[1,1],[0,mstrain],'k');
plot((1:c3),0*(1:c3),'k');

subplot(2,1,2)
plot(v0*visc0/charLength*vShear(1:c3)-eShear(1:c3))
str = ['Stress, \omega=', num2str(w)];
title(str);
hold on
[mstress,tstress]=max(v0*visc0/charLength*vShear-eShear);
plot(tstress*[1,1],[0,mstress],'k');
plot((1:c3),0*(1:c3),'k');

deltaw=(tstrain-tstress)*dt*w;
G1w=mstress/mstrain*cos(deltaw);
G2w=mstress/mstrain*sin(deltaw);

if(exist('text_handle_1') ==1)
    delete(text_handle_1);
end

str=['G_1^{\prime}=',num2str(G1w), ' G_2^{\prime\prime}=',num2str(G2w),' \delta=',num2str(deltaw)];
text_handle_1=text(10,mstress,str);
