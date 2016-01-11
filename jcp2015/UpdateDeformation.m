function [px,py,pz]=UpdateDeformation(px,py,pz,ux,uy,uz,dt,dx,charLength,v0)

[Em,En,Ep]=size(ux);

sx=zeros(size(ux));
sy=zeros(size(uy));
sz=zeros(size(uz));

sx(ux>0)=1;
sy(uy>0)=1;
sz(uz>0)=1;
sx(ux<0)=-1;
sy(uy<0)=-1;
sz(uz<0)=-1;

L=charLength;

px(2:Em-1,2:En-1,2:Ep-1)=px(2:Em-1,2:En-1,2:Ep-1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,2:En-1,2:Ep-1)+1).*ux(2:Em-1,2:En-1,2:Ep-1).*(3*px(2:Em-1,3:En,2:Ep-1)-4*px(2:Em-1,2:En-1,2:Ep-1)+px(2:Em-1,1:En-2,2:Ep-1)))+...
    ((sy(2:Em-1,2:En-1,2:Ep-1)+1).*uy(2:Em-1,2:En-1,2:Ep-1).*(3*px(3:Em,2:En-1,2:Ep-1)-4*px(2:Em-1,2:En-1,2:Ep-1)+px(1:Em-2,2:En-1,2:Ep-1)))+...
    ((sz(2:Em-1,2:En-1,2:Ep-1)+1).*uz(2:Em-1,2:En-1,2:Ep-1).*(3*px(2:Em-1,2:En-1,3:Ep)-4*px(2:Em-1,2:En-1,2:Ep-1)+px(2:Em-1,2:En-1,1:Ep-2)))+...
    ((sx(2:Em-1,2:En-1,2:Ep-1)-1).*ux(2:Em-1,2:En-1,2:Ep-1).*(px(2:Em-1,3:En,2:Ep-1)-4*px(2:Em-1,2:En-1,2:Ep-1)+3*px(2:Em-1,1:En-2,2:Ep-1)))+...
    ((sy(2:Em-1,2:En-1,2:Ep-1)-1).*uy(2:Em-1,2:En-1,2:Ep-1).*(px(3:Em,2:En-1,2:Ep-1)-4*px(2:Em-1,2:En-1,2:Ep-1)+3*px(1:Em-2,2:En-1,2:Ep-1)))+...
    ((sz(2:Em-1,2:En-1,2:Ep-1)-1).*uz(2:Em-1,2:En-1,2:Ep-1).*(px(2:Em-1,2:En-1,3:Ep)-4*px(2:Em-1,2:En-1,2:Ep-1)+3*px(2:Em-1,2:En-1,1:Ep-2))));

py(2:Em-1,2:En-1,2:Ep-1)=py(2:Em-1,2:En-1,2:Ep-1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,2:En-1,2:Ep-1)+1).*ux(2:Em-1,2:En-1,2:Ep-1).*(3*py(2:Em-1,3:En,2:Ep-1)-4*py(2:Em-1,2:En-1,2:Ep-1)+py(2:Em-1,1:En-2,2:Ep-1)))+...
    ((sy(2:Em-1,2:En-1,2:Ep-1)+1).*uy(2:Em-1,2:En-1,2:Ep-1).*(3*py(3:Em,2:En-1,2:Ep-1)-4*py(2:Em-1,2:En-1,2:Ep-1)+py(1:Em-2,2:En-1,2:Ep-1)))+...
    ((sz(2:Em-1,2:En-1,2:Ep-1)+1).*uz(2:Em-1,2:En-1,2:Ep-1).*(3*py(2:Em-1,2:En-1,3:Ep)-4*py(2:Em-1,2:En-1,2:Ep-1)+py(2:Em-1,2:En-1,1:Ep-2)))+...
    ((sx(2:Em-1,2:En-1,2:Ep-1)-1).*ux(2:Em-1,2:En-1,2:Ep-1).*(py(2:Em-1,3:En,2:Ep-1)-4*py(2:Em-1,2:En-1,2:Ep-1)+3*py(2:Em-1,1:En-2,2:Ep-1)))+...
    ((sy(2:Em-1,2:En-1,2:Ep-1)-1).*uy(2:Em-1,2:En-1,2:Ep-1).*(py(3:Em,2:En-1,2:Ep-1)-4*py(2:Em-1,2:En-1,2:Ep-1)+3*py(1:Em-2,2:En-1,2:Ep-1)))+...
    ((sz(2:Em-1,2:En-1,2:Ep-1)-1).*uz(2:Em-1,2:En-1,2:Ep-1).*(py(2:Em-1,2:En-1,3:Ep)-4*py(2:Em-1,2:En-1,2:Ep-1)+3*py(2:Em-1,2:En-1,1:Ep-2))));

pz(2:Em-1,2:En-1,2:Ep-1)=pz(2:Em-1,2:En-1,2:Ep-1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,2:En-1,2:Ep-1)+1).*ux(2:Em-1,2:En-1,2:Ep-1).*(3*pz(2:Em-1,3:En,2:Ep-1)-4*pz(2:Em-1,2:En-1,2:Ep-1)+pz(2:Em-1,1:En-2,2:Ep-1)))+...
    ((sy(2:Em-1,2:En-1,2:Ep-1)+1).*uy(2:Em-1,2:En-1,2:Ep-1).*(3*pz(3:Em,2:En-1,2:Ep-1)-4*pz(2:Em-1,2:En-1,2:Ep-1)+pz(1:Em-2,2:En-1,2:Ep-1)))+...
    ((sz(2:Em-1,2:En-1,2:Ep-1)+1).*uz(2:Em-1,2:En-1,2:Ep-1).*(3*pz(2:Em-1,2:En-1,3:Ep)-4*pz(2:Em-1,2:En-1,2:Ep-1)+pz(2:Em-1,2:En-1,1:Ep-2)))+...
    ((sx(2:Em-1,2:En-1,2:Ep-1)-1).*ux(2:Em-1,2:En-1,2:Ep-1).*(pz(2:Em-1,3:En,2:Ep-1)-4*pz(2:Em-1,2:En-1,2:Ep-1)+3*pz(2:Em-1,1:En-2,2:Ep-1)))+...
    ((sy(2:Em-1,2:En-1,2:Ep-1)-1).*uy(2:Em-1,2:En-1,2:Ep-1).*(pz(3:Em,2:En-1,2:Ep-1)-4*pz(2:Em-1,2:En-1,2:Ep-1)+3*pz(1:Em-2,2:En-1,2:Ep-1)))+...
    ((sz(2:Em-1,2:En-1,2:Ep-1)-1).*uz(2:Em-1,2:En-1,2:Ep-1).*(pz(2:Em-1,2:En-1,3:Ep)-4*pz(2:Em-1,2:En-1,2:Ep-1)+3*pz(2:Em-1,2:En-1,1:Ep-2))));




px(2:Em-1,1,2:Ep-1)=px(2:Em-1,1,2:Ep-1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,1,2:Ep-1)+1).*ux(2:Em-1,1,2:Ep-1).*(3*px(2:Em-1,2,2:Ep-1)-4*px(2:Em-1,1,2:Ep-1)+px(2:Em-1,En-1,2:Ep-1)))+...
    ((sy(2:Em-1,1,2:Ep-1)+1).*uy(2:Em-1,1,2:Ep-1).*(3*px(3:Em,1,2:Ep-1)-4*px(2:Em-1,1,2:Ep-1)+px(1:Em-2,1,2:Ep-1)))+...
    ((sz(2:Em-1,1,2:Ep-1)+1).*uz(2:Em-1,1,2:Ep-1).*(3*px(2:Em-1,1,3:Ep)-4*px(2:Em-1,1,2:Ep-1)+px(2:Em-1,1,1:Ep-2)))+...
    ((sx(2:Em-1,1,2:Ep-1)-1).*ux(2:Em-1,1,2:Ep-1).*(px(2:Em-1,2,2:Ep-1)-4*px(2:Em-1,1,2:Ep-1)+3*px(2:Em-1,En-1,2:Ep-1)))+...
    ((sy(2:Em-1,1,2:Ep-1)-1).*uy(2:Em-1,1,2:Ep-1).*(px(3:Em,1,2:Ep-1)-4*px(2:Em-1,1,2:Ep-1)+3*px(1:Em-2,1,2:Ep-1)))+...
    ((sz(2:Em-1,1,2:Ep-1)-1).*uz(2:Em-1,1,2:Ep-1).*(px(2:Em-1,1,3:Ep)-4*px(2:Em-1,1,2:Ep-1)+3*px(2:Em-1,1,1:Ep-2))));

py(2:Em-1,1,2:Ep-1)=py(2:Em-1,1,2:Ep-1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,1,2:Ep-1)+1).*ux(2:Em-1,1,2:Ep-1).*(3*py(2:Em-1,2,2:Ep-1)-4*py(2:Em-1,1,2:Ep-1)+py(2:Em-1,En-1,2:Ep-1)))+...
    ((sy(2:Em-1,1,2:Ep-1)+1).*uy(2:Em-1,1,2:Ep-1).*(3*py(3:Em,1,2:Ep-1)-4*py(2:Em-1,1,2:Ep-1)+py(1:Em-2,1,2:Ep-1)))+...
    ((sz(2:Em-1,1,2:Ep-1)+1).*uz(2:Em-1,1,2:Ep-1).*(3*py(2:Em-1,1,3:Ep)-4*py(2:Em-1,1,2:Ep-1)+py(2:Em-1,1,1:Ep-2)))+...
    ((sx(2:Em-1,1,2:Ep-1)-1).*ux(2:Em-1,1,2:Ep-1).*(py(2:Em-1,2,2:Ep-1)-4*py(2:Em-1,1,2:Ep-1)+3*py(2:Em-1,En-1,2:Ep-1)))+...
    ((sy(2:Em-1,1,2:Ep-1)-1).*uy(2:Em-1,1,2:Ep-1).*(py(3:Em,1,2:Ep-1)-4*py(2:Em-1,1,2:Ep-1)+3*py(1:Em-2,1,2:Ep-1)))+...
    ((sz(2:Em-1,1,2:Ep-1)-1).*uz(2:Em-1,1,2:Ep-1).*(py(2:Em-1,1,3:Ep)-4*py(2:Em-1,1,2:Ep-1)+3*py(2:Em-1,1,1:Ep-2))));

pz(2:Em-1,1,2:Ep-1)=pz(2:Em-1,1,2:Ep-1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,1,2:Ep-1)+1).*ux(2:Em-1,1,2:Ep-1).*(3*pz(2:Em-1,2,2:Ep-1)-4*pz(2:Em-1,1,2:Ep-1)+pz(2:Em-1,En-1,2:Ep-1)))+...
    ((sy(2:Em-1,1,2:Ep-1)+1).*uy(2:Em-1,1,2:Ep-1).*(3*pz(3:Em,1,2:Ep-1)-4*pz(2:Em-1,1,2:Ep-1)+pz(1:Em-2,1,2:Ep-1)))+...
    ((sz(2:Em-1,1,2:Ep-1)+1).*uz(2:Em-1,1,2:Ep-1).*(3*pz(2:Em-1,1,3:Ep)-4*pz(2:Em-1,1,2:Ep-1)+pz(2:Em-1,1,1:Ep-2)))+...
    ((sx(2:Em-1,1,2:Ep-1)-1).*ux(2:Em-1,1,2:Ep-1).*(pz(2:Em-1,2,2:Ep-1)-4*pz(2:Em-1,1,2:Ep-1)+3*pz(2:Em-1,En-1,2:Ep-1)))+...
    ((sy(2:Em-1,1,2:Ep-1)-1).*uy(2:Em-1,1,2:Ep-1).*(pz(3:Em,1,2:Ep-1)-4*pz(2:Em-1,1,2:Ep-1)+3*pz(1:Em-2,1,2:Ep-1)))+...
    ((sz(2:Em-1,1,2:Ep-1)-1).*uz(2:Em-1,1,2:Ep-1).*(pz(2:Em-1,1,3:Ep)-4*pz(2:Em-1,1,2:Ep-1)+3*pz(2:Em-1,1,1:Ep-2))));

px(:,En,:)=px(:,1,:);
py(:,En,:)=py(:,1,:);
pz(:,En,:)=pz(:,1,:);


px(2:Em-1,2:En-1,1)=px(2:Em-1,2:En-1,1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,2:En-1,1)+1).*ux(2:Em-1,2:En-1,1).*(3*px(2:Em-1,3:En,1)-4*px(2:Em-1,2:En-1,1)+px(2:Em-1,1:En-2,1)))+...
    ((sy(2:Em-1,2:En-1,1)+1).*uy(2:Em-1,2:En-1,1).*(3*px(3:Em,2:En-1,1)-4*px(2:Em-1,2:En-1,1)+px(1:Em-2,2:En-1,1)))+...
    ((sz(2:Em-1,2:En-1,1)+1).*uz(2:Em-1,2:En-1,1).*(3*px(2:Em-1,2:En-1,2)-4*px(2:Em-1,2:En-1,1)+px(2:Em-1,2:En-1,Ep-1)))+...
    ((sx(2:Em-1,2:En-1,1)-1).*ux(2:Em-1,2:En-1,1).*(px(2:Em-1,3:En,1)-4*px(2:Em-1,2:En-1,1)+3*px(2:Em-1,1:En-2,1)))+...
    ((sy(2:Em-1,2:En-1,1)-1).*uy(2:Em-1,2:En-1,1).*(px(3:Em,2:En-1,1)-4*px(2:Em-1,2:En-1,1)+3*px(1:Em-2,2:En-1,1)))+...
    ((sz(2:Em-1,2:En-1,1)-1).*uz(2:Em-1,2:En-1,1).*(px(2:Em-1,2:En-1,2)-4*px(2:Em-1,2:En-1,1)+3*px(2:Em-1,2:En-1,Ep-1))));

py(2:Em-1,2:En-1,1)=py(2:Em-1,2:En-1,1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,2:En-1,1)+1).*ux(2:Em-1,2:En-1,1).*(3*py(2:Em-1,3:En,1)-4*py(2:Em-1,2:En-1,1)+py(2:Em-1,1:En-2,1)))+...
    ((sy(2:Em-1,2:En-1,1)+1).*uy(2:Em-1,2:En-1,1).*(3*py(3:Em,2:En-1,1)-4*py(2:Em-1,2:En-1,1)+py(1:Em-2,2:En-1,1)))+...
    ((sz(2:Em-1,2:En-1,1)+1).*uz(2:Em-1,2:En-1,1).*(3*py(2:Em-1,2:En-1,2)-4*py(2:Em-1,2:En-1,1)+py(2:Em-1,2:En-1,Ep-1)))+...
    ((sx(2:Em-1,2:En-1,1)-1).*ux(2:Em-1,2:En-1,1).*(py(2:Em-1,3:En,1)-4*py(2:Em-1,2:En-1,1)+3*py(2:Em-1,1:En-2,1)))+...
    ((sy(2:Em-1,2:En-1,1)-1).*uy(2:Em-1,2:En-1,1).*(py(3:Em,2:En-1,1)-4*py(2:Em-1,2:En-1,1)+3*py(1:Em-2,2:En-1,1)))+...
    ((sz(2:Em-1,2:En-1,1)-1).*uz(2:Em-1,2:En-1,1).*(py(2:Em-1,2:En-1,2)-4*py(2:Em-1,2:En-1,1)+3*py(2:Em-1,2:En-1,Ep-1))));

pz(2:Em-1,2:En-1,1)=pz(2:Em-1,2:En-1,1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,2:En-1,1)+1).*ux(2:Em-1,2:En-1,1).*(3*pz(2:Em-1,3:En,1)-4*pz(2:Em-1,2:En-1,1)+pz(2:Em-1,1:En-2,1)))+...
    ((sy(2:Em-1,2:En-1,1)+1).*uy(2:Em-1,2:En-1,1).*(3*pz(3:Em,2:En-1,1)-4*pz(2:Em-1,2:En-1,1)+pz(1:Em-2,2:En-1,1)))+...
    ((sz(2:Em-1,2:En-1,1)+1).*uz(2:Em-1,2:En-1,1).*(3*pz(2:Em-1,2:En-1,2)-4*pz(2:Em-1,2:En-1,1)+pz(2:Em-1,2:En-1,Ep-1)))+...
    ((sx(2:Em-1,2:En-1,1)-1).*ux(2:Em-1,2:En-1,1).*(pz(2:Em-1,3:En,1)-4*pz(2:Em-1,2:En-1,1)+3*pz(2:Em-1,1:En-2,1)))+...
    ((sy(2:Em-1,2:En-1,1)-1).*uy(2:Em-1,2:En-1,1).*(pz(3:Em,2:En-1,1)-4*pz(2:Em-1,2:En-1,1)+3*pz(1:Em-2,2:En-1,1)))+...
    ((sz(2:Em-1,2:En-1,1)-1).*uz(2:Em-1,2:En-1,1).*(pz(2:Em-1,2:En-1,2)-4*pz(2:Em-1,2:En-1,1)+3*pz(2:Em-1,2:En-1,Ep-1))));

px(:,:,Ep)=px(:,:,1);
py(:,:,Ep)=py(:,:,1);
pz(:,:,Ep)=pz(:,:,1);



px(2:Em-1,1,1)=px(2:Em-1,1,1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,1,1)+1).*ux(2:Em-1,1,1).*(3*px(2:Em-1,2,1)-4*px(2:Em-1,1,1)+px(2:Em-1,En-1,1)))+...
    ((sy(2:Em-1,1,1)+1).*uy(2:Em-1,1,1).*(3*px(3:Em,1,1)-4*px(2:Em-1,1,1)+px(1:Em-2,1,1)))+...
    ((sz(2:Em-1,1,1)+1).*uz(2:Em-1,1,1).*(3*px(2:Em-1,1,2)-4*px(2:Em-1,1,1)+px(2:Em-1,1,Ep-1)))+...
    ((sx(2:Em-1,1,1)-1).*ux(2:Em-1,1,1).*(px(2:Em-1,2,1)-4*px(2:Em-1,1,1)+3*px(2:Em-1,En-1,1)))+...
    ((sy(2:Em-1,1,1)-1).*uy(2:Em-1,1,1).*(px(3:Em,1,1)-4*px(2:Em-1,1,1)+3*px(1:Em-2,1,1)))+...
    ((sz(2:Em-1,1,1)-1).*uz(2:Em-1,1,1).*(px(2:Em-1,1,2)-4*px(2:Em-1,1,1)+3*px(2:Em-1,1,Ep-1))));

py(2:Em-1,1,1)=py(2:Em-1,1,1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,1,1)+1).*ux(2:Em-1,1,1).*(3*py(2:Em-1,2,1)-4*py(2:Em-1,1,1)+py(2:Em-1,En-1,1)))+...
    ((sy(2:Em-1,1,1)+1).*uy(2:Em-1,1,1).*(3*py(3:Em,1,1)-4*py(2:Em-1,1,1)+py(1:Em-2,1,1)))+...
    ((sz(2:Em-1,1,1)+1).*uz(2:Em-1,1,1).*(3*py(2:Em-1,1,2)-4*py(2:Em-1,1,1)+py(2:Em-1,1,Ep-1)))+...
    ((sx(2:Em-1,1,1)-1).*ux(2:Em-1,1,1).*(py(2:Em-1,2,1)-4*py(2:Em-1,1,1)+3*py(2:Em-1,En-1,1)))+...
    ((sy(2:Em-1,1,1)-1).*uy(2:Em-1,1,1).*(py(3:Em,1,1)-4*py(2:Em-1,1,1)+3*py(1:Em-2,1,1)))+...
    ((sz(2:Em-1,1,1)-1).*uz(2:Em-1,1,1).*(py(2:Em-1,1,2)-4*py(2:Em-1,1,1)+3*py(2:Em-1,1,Ep-1))));

pz(2:Em-1,1,1)=pz(2:Em-1,1,1)+0.5*0.5*dt/dx*v0/L*(...
    ((sx(2:Em-1,1,1)+1).*ux(2:Em-1,1,1).*(3*pz(2:Em-1,2,1)-4*pz(2:Em-1,1,1)+pz(2:Em-1,En-1,1)))+...
    ((sy(2:Em-1,1,1)+1).*uy(2:Em-1,1,1).*(3*pz(3:Em,1,1)-4*pz(2:Em-1,1,1)+pz(1:Em-2,1,1)))+...
    ((sz(2:Em-1,1,1)+1).*uz(2:Em-1,1,1).*(3*pz(2:Em-1,1,2)-4*pz(2:Em-1,1,1)+pz(2:Em-1,1,Ep-1)))+...
    ((sx(2:Em-1,1,1)-1).*ux(2:Em-1,1,1).*(pz(2:Em-1,2,1)-4*pz(2:Em-1,1,1)+3*pz(2:Em-1,En-1,1)))+...
    ((sy(2:Em-1,1,1)-1).*uy(2:Em-1,1,1).*(pz(3:Em,1,1)-4*pz(2:Em-1,1,1)+3*pz(1:Em-2,1,1)))+...
    ((sz(2:Em-1,1,1)-1).*uz(2:Em-1,1,1).*(pz(2:Em-1,1,2)-4*pz(2:Em-1,1,1)+3*pz(2:Em-1,1,Ep-1))));

px(:,En,Ep)=px(:,1,1);
px(:,En,1)=px(:,1,1);
px(:,1,Ep)=px(:,1,1);

py(:,En,Ep)=py(:,1,1);
py(:,En,1)=py(:,1,1);
py(:,1,Ep)=py(:,1,1);

pz(:,En,Ep)=pz(:,1,1);
pz(:,En,1)=pz(:,1,1);
pz(:,1,Ep)=pz(:,1,1);

%pz(Em,:,:)=pz(Em,:,:)=0.5*0.5*dt/dx*v0/L*uz(Em,:,:)*(

end