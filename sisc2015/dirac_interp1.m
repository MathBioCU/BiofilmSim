function val=dirac_interp1(r)
%every value in r should be between 0 and 2
val=r;
rless=r(r<=1);
rgreat=r(r>1 & r<=2);
% val=zeros(length(r),1);
val(r<=1)=(3-2*rless+sqrt(1+4*rless-4*rless.^2))/8;
val(r>1 & r<=2)=(5-2*rgreat-sqrt(-7+12*rgreat-4*rgreat.^2))/8;
val(r>2)=0;
% val(absr>2)=0;