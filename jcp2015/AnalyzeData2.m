[m,t]=max(fStrain);
[m2,t2]=max(v0*visc0/charLength*vShear-eShear);
delta=(t-t2)*dt*w;
G1=m2/m*cos(delta);
G2=m2/m*sin(delta);

fprintf('G1 %4.4f  G2 %4.4f  delta %4.4f\n', G1,G2,delta);