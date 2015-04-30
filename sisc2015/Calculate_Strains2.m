function [bStrain, fStrain,S] = Calculate_Strains2(A,X,S,ux,uy,uz,dx,X0,S0,dt,st,xlength,ylength,zlength,upper,h,x,y,z,numOfnonzero, bthickness)

%[Em,En,Ep]=size(ux);
%Emm=Em-1;
%Emmm=Em-2;
%Epm=Ep-1;
%Epmm=Ep-2;
%n=length(X(:,1));


%other idea for getting biofilm strain- interpolate U to Eulerian grid
%using transfer function, and then calculate the same way the fluid strain
%is being calculated.s
U=X-X0; %displacement
Xr=[X(:,1),X(:,2),X(:,3)+zlength];
Xl=[X(:,1),X(:,2),X(:,3)-zlength];
Ur=Xr-X0;
Ul=Xl-X0;
Xr=[X(:,1)+xlength,X(:,2),X(:,3)];
Xl=[X(:,1)-xlength,X(:,2),X(:,3)];
Uf=Xr-X0;
Ub=Xl-X0;

%zind=find(abs(U)>abs(Ur));
%U(zind)=Ur(zind);
%zind=find(abs(U)>abs(Ul));
%U(zind)=Ul(zind);
%xind=find(abs(U)>abs(Uf));
%U(xind)=Uf(xind);
%xind=find(abs(U)>abs(Ub));
%U(xind)=Ub(xind);

%EUz=transferLtoEDisp3Dper_e2(h,x,y,z,X,U(:,3),numOfnonzero,zlength,xlength,1/20);
%EUz=EUz*max(max(max(abs(U(:,3)))))/max(max(max(abs(EUz))));
%EUy=transferLtoEDisp3Dper_e2(h,x,y,z,X,U(:,2),numOfnonzero,zlength,xlength,1/20);
%EUy=EUy*max(max(max(abs(U(:,2)))))/max(max(max(abs(EUy))));

% EUx=transferLtoEDisp3Dper_e(h,x,y,z,X,U(:,1),numOfnonzero,zlength);
% EUx=EUx*max(max(max(abs(U(:,1)))))/max(max(max(abs(EUx))));


%Shear_Strainb=1/2*((EUz(3:Em,2:En,2:Epm)-EUz(1:Emmm,2:En,2:Epm))/(2*dx)+(EUy(2:Emm,2:En,3:Ep)-EUy(2:Emm,2:En,1:Epmm))/(2*dx));


%bStrain=1/2*((sum(sum(EUz(Em,2:En,1:Epm)-EUz(Emm,2:En,1:Epm)))/(dx)+sum(sum(-3*EUz(1,2:En,1:Epm)+4*EUz(2,2:En,1:Epm)-EUz(3,2:En,1:Epm)))/(2*dx)...
%    +sum(sum(sum(EUz(3:Em,2:En,1:Epm)-EUz(1:Emmm,2:En,1:Epm))))/(2*dx))*dx^3/(xlength*ylength*zlength)+sum(sum(sum((EUy(:,2:En,3:Ep)-EUy(:,2:En,1:Epmm))/(2*dx))))*dx^3/(xlength*ylength*zlength));

%bStraintop=1/2*((...
%    +sum(sum(sum(EUz(Em-4:Em-2,2:En,1:Epm)-EUz(Em-6:Em-4,2:En,1:Epm))))/(2*dx))*dx^3/(xlength*zlength*h*2)+sum(sum(sum((EUy(Em-4:Em-2,2:En,3:Ep)-EUy(Em-6:Em-4,2:En,1:Epmm))/(2*dx))))*dx^3/(xlength*2*h*zlength));



strain=zeros(length(upper),1);
for i=1:length(upper)
	ind=find(A(upper(i),:)~=0);
	ind2=find(X(ind,2)<ylength-3*h);
	if length(ind2)>0
		for j=1:length(ind2)
			strain(i)=strain(i)+1/length(ind2)*(U(upper(i),3)-U(ind(ind2(j)),3))/(X(upper(i),2)-X(ind(ind2(j)),2));
		end
	end
end




bStrain = 1/2*mean(strain);

S1=S;
S1(:,3)=mod(S1(:,3),zlength);
S1(:,1)=mod(S1(:,1),xlength);
S1(S1(:,2)>ylength,2)=ylength;

Us=transferEtoLvel3Dper_e2a(h,ux,uy,uz,x,y,z,S1,zlength,xlength);
S=S+dt/st*Us;

V=S-S0; %displacement

Sr=[S(:,1),S(:,2),S(:,3)+zlength];
Sl=[S(:,1),S(:,2),S(:,3)-zlength];
Vr=Sr-S0;
Vl=Sl-S0;
Sr=[S(:,1)+xlength,S(:,2),S(:,3)];
Sl=[S(:,1)-xlength,S(:,2),S(:,3)];
Vf=Sr-S0;
Vb=Sl-S0;


%zind=find(abs(V)>abs(Vr));
%V(zind)=Vr(zind);
%zind=find(abs(V)>abs(Vl));
%V(zind)=Vl(zind);
%xind=find(abs(V)>abs(Vf));
%V(xind)=Vf(xind);
%xind=find(abs(V)>abs(Vb));
%V(xind)=Vb(xind);

LVz=zeros(100,1);
for i=1:100
    h1=S(i,2)-S(i+100,2);
    h2=S(i+100,2)-S(i+200,2);
%     LVz(i)=1/2*(h1+h2)*(h2*(V(i,3)-V(i+100,3))/(h1)+h1*(V(i+100,3)-V(i+200,3))/(h2));
    LVz(i)=1/2*((V(i,3)-V(i+100,3))/(h1)+(V(i+100,3)-V(i+200,3))/(h2));
end
% ind=find(LVz<mean(LVz)+std(LVz) & LVz>mean(LVz)-std(LVz));
% LVz=LVz(ind);
%EVz=transferLtoEDisp3Dper_e2(h,x,y,z,S(37:72,:),LVz,36,zlength,xlength,2*h);
%EVz=EVz*max(max(max(abs(LVz))))/max(max(max(abs(EVz))));
%EVy=transferLtoEDisp3Dper_e2(h,x,y,z,S,V(:,2),36,zlength,xlength,2*h);
%EVy=EVy*max(max(max(abs(V(:,2)))))/max(max(max(abs(EVy))));

%fStraintop=sum(sum(sum(EVz)))*dx^3*(xlength/h*zlength/h)/36*1/(xlength*zlength*h);
%fStraintop=1;
fStrain=1/2*mean(LVz);
%Shear_Strainf=1;


%%select points near top. ind=find(..) gives the 1-dimensional index of these points in the 3D array.
%ind=find(y>=1.8-2*h);
%Xt=zeros(length(ind),3);
%Xt(:,1)=x(ind);
%Xt(:,2)=y(ind);
%Xt(:,3)=z(ind);
%%select Lagrangian positions that will be included in the compact support of Delta functions around indexed points above
%ind=find(S{2}>=1.8-3*h);

%St=zeros(length(ind),3);
%St(:,1)=mod(S{1}(ind),xlength);
%St(:,2)=S{2}(ind2);
%St(:,3)=mod(S{3}(ind),zlength);

%toc
%find labels for points
%Ax=transferLtoEDisp3Dper_e2(h,x,y,z,St,Xt(:,1),length(Xt),zlength,xlength,2*h);
%Ax=Ax*max(max(max(abs(Ax))))/max(max(max(abs(S{1}))));

%Ay=transferLtoEDisp3Dper_e2(h,x,y,z,St,Xt(:,2),length(Xt),zlength,xlength,2*h);
%Ay=Ay*max(max(max(abs(Ay))))/max(max(max(abs(S{2}))));
%Az=transferLtoEDisp3Dper_e2(h,x,y,z,St,Xt(:,3),length(Xt),zlength,xlength,2*h);
%Az=Az*max(max(max(abs(Az))))/max(max(max(abs(S{3}))));
%toc
%convert labels to 3D array (currently inefficient)

%fStrain=1/2*(sum(sum(3*S{3}(Em,2:En,1:Epm)-4*S{3}(Emm,2:En,1:Epm)+S{3}(Emmm,2:En,1:Epm)))/(2*dx)+sum(sum(-3*S{3}(1,2:En,1:Epm)+4*S{3}(2,2:En,1:Epm)-S{3}(3,2:En,1:Epm)))/(2*dx)...
%    +sum(sum(sum(S{3}(3:Em,2:En,1:Epm)-S{3}(1:Emmm,2:En,1:Epm))))/(2*dx))*dx^3/(xlength*ylength*zlength);

%fStraintop=1/2*(sum(sum(3*S{3}(Em,2:En,1:Epm)-4*S{3}(Emm,2:En,1:Epm)+S{3}(Emmm,2:En,1:Epm)))/(2*dx)...
%    +sum(sum(sum(S{3}(Em-5:Em,2:En,1:Epm)-S{3}(Em-7:Emmm,2:En,1:Epm))))/(2*dx))*dx^3/(xlength*ylength*h*5);

%fStraintop=1/2*(sum(sum((2*dx)/(3*Az(3,2:En,1:Epm)-4*Az(2,2:En,1:Epm)+Az(1,2:En,1:Epm))))...
%    +sum(sum(sum((2*dx)/(Az(1:2,2:En,1:Epm)-Az(3:4,2:En,1:Epm))))))*dx^3/(xlength*ylength*h*3);


%compute strain at top

%EU{1}=EUz;
%EU{2}=EUy;
% EU{3}=EUx;


%Shear_Strainf=1/2*((S{3}(3:Em,2:En,1:Epm)-S{3}(1:Emmm,2:En,1:Epm))/(2*dx));

% bStrain33=(EUz(:,:,3:Ep)-EUz(:,:,1:Ep-2))/(2*dx);
% bStrain32=1/2*((EUz(3:Em,:,2:Epm)-EUz(1:Em-2,:,2:Epm))/(2*dx)+(EUy(2:Emm,:,3:Ep)-EUy(2:Emm,:,1:Ep-2))/(2*dx));
% bStrain31=1/2*((EUz(:,3:En,2:Epm)-EUz(:,1:En-2,2:Epm))/(2*dx)+(EUx(:,2:En-1,3:Ep)-EUx(:,2:En-1,1:Ep-2))/(2*dx));
% bStrain22=(EUy(3:Em,:,:)-EUy(1:Em-2,:,:))/(2*dx);
% bStrain21=1/2*((EUy(2:Emm,3:En,:)-EUy(2:Emm,1:En-2,:))/(2*dx)+(EUx(3:Em,2:En-1,:)-EUx(1:Em-2,2:En-1,:))/(2*dx));
% bStrain11=(EUx(:,3:En,:)-EUx(:,1:En-2,:))/(2*dx);
% 
% tbStrain=cell(6,1);
% tbStrain{1}=bStrain11; tbStrain{2}=bStrain22; tbStrain{3}=bStrain33;
% tbStrain{4}=bStrain21; tbStrain{5}=bStrain32; tbStrain{6}=bStrain31;
% 
% fStrain33=(S{3}(:,:,3:Ep)-S{3}(:,:,1:Ep-2))/(2*dx);
% fStrain32=1/2*((S{3}(3:Em,:,2:Epm)-S{3}(1:Em-2,:,2:Epm))/(2*dx)+(S{2}(2:Emm,:,3:Ep)-S{2}(2:Emm,:,1:Ep-2))/(2*dx));
% fStrain31=1/2*((S{3}(:,3:En,2:Epm)-S{3}(:,1:En-2,2:Epm))/(2*dx)+(S{1}(:,2:En-1,3:Ep)-S{1}(:,2:En-1,1:Ep-2))/(2*dx));
% fStrain22=(S{2}(3:Em,:,:)-S{2}(1:Em-2,:,:))/(2*dx);
% fStrain21=1/2*((S{2}(2:Emm,3:En,:)-S{2}(2:Emm,1:En-2,:))/(2*dx)+(S{1}(3:Em,2:En-1,:)-S{1}(1:Em-2,2:En-1,:))/(2*dx));
% fStrain11=(S{1}(:,3:En,:)-S{1}(:,1:En-2,:))/(2*dx);
% 
% tfStrain=cell(6,1);
% tfStrain{1}=fStrain11; tfStrain{2}=fStrain22; tfStrain{3}=fStrain33;
% tfStrain{4}=fStrain21; tfStrain{5}=fStrain32; tfStrain{6}=fStrain31;


end

