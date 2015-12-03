clc
%clear all 

%profile on
w=4.991;
B01=0;

uxx=cell(1,4);
uyy=uxx;
uzz=uxx;
Xx=uxx;
Pp=uxx;

for ii=1:4
    %w=fr(i);
    b=B01;
    E=0.000000;

    dt=1/(250*2^(ii-1));
    dt=dt/w;
    numtimesteps=125*2^(ii-1);
    charLength=10*10^-6;
    fmax=25000;
    addlvisc=0;
    connectdist=3/18.5;
    sigma0=0;
    JMain3DSim7
    uxx{ii}=ux;
    uyy{ii}=uy;
    uzz{ii}=uz;
    Pp{ii}=p;
    Xx{ii}=X;

    [G1,G2,Delta,max_stress,min_stress]=Analyze_Data2_3D(w,w*dt,e0,v0*visc0/charLength*vShear-eShear,fStrain,numtimesteps);
    str1=num2str(w);
    str2=num2str(fmax/1000);
    str3=num2str(B01);
    str4=num2str(1/(w*dt));
    str5=num2str(connectdist);
    str6=num2str(addlvisc);
    str7=num2str(sigma0);
    runid=['w',str1,'_f',str2,'_b',str3,'_cnd',str5,'_visc',str6,'_dt',str4,'_conv'];
    save([runid,'.mat'])
    Complete='Complete';
    G1s=num2str(G1);
    G2s=num2str(G2);
    maxstrain=num2str(max(fStraintop+bStraintop));
    fprintf('G1 %s,G2  %s,Max Strain %s, %s\n',G1s,G2s,maxstrain,Complete);	

end

 uz_exact=zeros(size(uz));
 kw=(w/(2*visc0/rho0))^(1/2);
 for i=1:Em
     uz_exact(i,:,:)=sqrt(real(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))).^2+imag(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))).^2).*sin(w*(dt*c3)+angle(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))));
 end


Eu1=sqrt(sum(sum(sum(((uxx{1}-uxx{2}).^2+(uyy{1}-uyy{2}).^2+(uzz{1}-uzz{2}).^2)*h^3))));
Eu2=sqrt(sum(sum(sum(((uxx{2}-uxx{3}).^2+(uyy{2}-uyy{3}).^2+(uzz{2}-uzz{3}).^2)*h^3))));
Eu3=sqrt(sum(sum(sum(((uxx{3}-uxx{4}).^2+(uyy{3}-uyy{4}).^2+(uzz{3}-uzz{4}).^2)*h^3))));

Ep1=sqrt(sum(sum(sum(((Pp{1}-Pp{2}).^2)*h^3))));
Ep2=sqrt(sum(sum(sum(((Pp{2}-Pp{3}).^2)*h^3))));
Ep3=sqrt(sum(sum(sum(((Pp{3}-Pp{4}).^2)*h^3))));

Ex1=sqrt(sum(sum(sum(((Xx{1}-Xx{2}).^2)*h^3))));
Ex2=sqrt(sum(sum(sum(((Xx{2}-Xx{3}).^2)*h^3))));
Ex3=sqrt(sum(sum(sum(((Xx{3}-Xx{4}).^2)*h^3))));

Ru1=log2(Eu1/Eu2);
Ru2=log2(Eu2/Eu3);
Ru=mean([Ru1,Ru2]);

Rp1=log2(Ep1/Ep2);
Rp2=log2(Ep2/Ep3);
Rp=mean([Rp1, Rp2]);

Rx1=log2(Ex1/Ex2);
Rx2=log2(Ex2/Ex3);
Rx=mean([Rx1,Rx2]);

%p=profile('info');
%save profile_results p;
%profile off
runid

