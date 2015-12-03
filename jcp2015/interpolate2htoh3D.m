function vh=interpolate2htoh3D(v2h)
[pm,pn,pp]=size(v2h);
Em=pm*2-1;
En=pn*2-1;
Ep=pp*2-1;
Emm=Em-1;
Enm=En-1;
Epm=Ep-1;
pmm=pm-1;
pnm=pn-1;
ppm=pp-1;
% p=log2(size1(1)-1);
vh=zeros(Em,En,Ep);
% h=2^(-p);
vh(1:2:Em,1:2:En,1:2:Ep)=v2h;
% v2h=[v2h,v2h(:,2)];

vh(1:2:Em,2:2:Enm,1:2:Ep)=(v2h(:,2:pn,:)+v2h(:,1:pnm,:))/2;
vh(1:2:Em,1:2:En,2:2:Epm)=(v2h(:,:,2:pp)+v2h(:,:,1:ppm))/2;
vh(2:2:Emm,1:2:En,1:2:Ep)=(v2h(2:pm,:,:)+v2h(1:pmm,:,:))/2;

vh(2:2:Emm,2:2:Enm,1:2:Ep)=(v2h(1:pmm,1:pnm,:)+v2h(1:pmm,2:pn,:)+v2h(2:pm,1:pnm,:)+v2h(2:pm,2:pn,:))/4;
vh(1:2:Em,2:2:Enm,2:2:Epm)=(v2h(:,1:pnm,1:ppm)+v2h(:,2:pn,1:ppm)+v2h(:,1:pnm,2:pp)+v2h(:,2:pn,2:pp))/4;
vh(2:2:Emm,1:2:En,2:2:Epm)=(v2h(1:pmm,:,1:ppm)+v2h(2:pm,:,1:ppm)+v2h(1:pmm,:,2:pp)+v2h(2:pm,:,2:pp))/4;

vh(2:2:Emm,2:2:Enm,2:2:Epm)=(v2h(1:pmm,1:pnm,1:ppm)+v2h(2:pm,1:pnm,1:ppm)+v2h(1:pmm,1:pnm,2:pp)+...
    v2h(2:pm,1:pnm,2:pp)+v2h(1:pmm,2:pn,1:ppm)+v2h(2:pm,2:pn,1:ppm)+v2h(1:pmm,2:pn,2:pp)+...
    v2h(2:pm,2:pn,2:pp))/8;


