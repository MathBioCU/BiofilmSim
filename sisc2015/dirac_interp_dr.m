function val=dirac_interp_dr(r)
%every value in r should be between 0 and 2
val=r;
rless1=r(r<0 & r>=-1);
rgreat1=r(r<-1 & r>=-2);
rless2=r(r>=0 & r<=1);
rgreat2=r(r>1 & r<=2);
% val=zeros(length(r),1);
val(r>=-1 & r<0)=(2+(-4-8*rless1)./(2*sqrt(1-4*rless1-4*rless1.^2)))/8;
val(r<-1 & r>=-2)=(2-(-12-8*rgreat1)./(2*sqrt(-7-12*rgreat1-4*rgreat1.^2)))/8;
val(r<=1 & r>=0)=(-2+(4-8*rless2)./(2*sqrt(1+4*rless2-4*rless2.^2)))/8;
val(r>1 & r<=2)=(-2-(12-8*rgreat2)./(2*sqrt(-7+12*rgreat2-4*rgreat2.^2)))/8;
val(abs(r)>2)=0;
% val(absr>2)=0;