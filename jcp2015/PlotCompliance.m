subplot(2,1,1)
loglog(times(1:c3),fStrain(1:c3))
str = ['Strain, \omega=', num2str(w)];
title(str);

subplot(2,1,2)
plot(times(1:c3),v0*visc0/charLength*vShear(1:c3)-eShear(1:c3))
str = ['Stress, \omega=', num2str(w)];
title(str);
