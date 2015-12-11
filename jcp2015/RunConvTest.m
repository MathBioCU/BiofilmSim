clc
clear

%profile on
w=49.91;
B01=0;

uxx=cell(1,4);
uyy=uxx;
uzz=uxx;
Xx=uxx;
Pp=uxx;

T=1;
S=0;

TS=[16, 47, 140,419];
if T==1
%temporal convergence tests
    for ii=1:4
        %w=fr(i);
        b=B01;
        E=0.000000;
        
        dx=1/32;
        dt=1/(250*3^(ii-1));
        dt=dt/w;
        numtimesteps=125*2^(ii-1);%TS(ii);
        charLength=10*10^-6;
        fmax=25000;
        addlvisc=500;
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
     %   save([runid,'.mat'])
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
    
    d0=sqrt(sum(sum(sum((uzz{1}-uz_exact).^2))));
    d1=sqrt(sum(sum(sum((uzz{1}-uz_exact).^2))));
    d2=sqrt(sum(sum(sum((uzz{1}-uz_exact).^2))));
    d3=sqrt(sum(sum(sum((uzz{1}-uz_exact).^2))));
    
    D1=log2(d0/d1);
    D2=log2(d1/d2);
    D3=log2(d2/d3);
    
    save('ConvergenceTestTimeBiofilm.mat');
end    
if S==1
    uz_exact=cell(1,4);

%spatial convergence test
    for ii=1:3
        %w=fr(i);
        b=B01;
        E=0.000000;
        dx=1/(16*2^(ii-1));
        dt=1/(500);
        dt=dt/w;
        numtimesteps=100;
        charLength=10*10^-6;
        fmax=25000;
        addlvisc=500;
        connectdist=3/18.5;
        sigma0=0;
        JMain3DSim7
        uxx{ii}=ux;
        uyy{ii}=uy;
        uzz{ii}=uz;
        Pp{ii}=p;
        Xx{ii}=X;
        
        uz_exact{ii}=zeros(size(uz));
        kw=(w/(2*visc0/rho0))^(1/2);
        for i=1:Em
            uz_exact{ii}(i,:,:)=sqrt(real(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))).^2+imag(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))).^2).*sin(w*(dt*c3)+angle(sinh(kw*y(i,:,:)*charLength*(1+sqrt(-1)))/sinh(kw*ylength*charLength*(1+sqrt(-1)))));
        end

        [G1,G2,Delta,max_stress,min_stress]=Analyze_Data2_3D(w,w*dt,e0,v0*visc0/charLength*vShear-eShear,fStrain,numtimesteps);
        str1=num2str(w);
        str2=num2str(fmax/1000);
        str3=num2str(B01);
        str4=num2str(1/(w*dt));
        str5=num2str(connectdist);
        str6=num2str(addlvisc);
        str7=num2str(sigma0);
        runid=['w',str1,'_f',str2,'_b',str3,'_cnd',str5,'_visc',str6,'_dt',str4,'_conv'];
      %  save([runid,'.mat'])
        Complete='Complete';
        G1s=num2str(G1);
        G2s=num2str(G2);
        maxstrain=num2str(max(fStraintop+bStraintop));
        fprintf('G1 %s,G2  %s,Max Strain %s, %s\n',G1s,G2s,maxstrain,Complete);	

    end


    Eu1=sqrt(sum(sum(sum(((interpolate2htoh3D(uxx{1})-uxx{2}).^2+(interpolate2htoh3D(uyy{1})-uyy{2}).^2+(interpolate2htoh3D(uzz{1})-uzz{2}).^2)*(4*h)^3))));
    Eu2=sqrt(sum(sum(sum(((interpolate2htoh3D(uxx{2})-uxx{3}).^2+(interpolate2htoh3D(uyy{2})-uyy{3}).^2+(interpolate2htoh3D(uzz{2})-uzz{3}).^2)*(2*h)^3))));
  %  Eu3=sqrt(sum(sum(sum(((interpolate2htoh3D(uxx{3})-uxx{4}).^2+(interpolate2htoh3D(uyy{3})-uyy{4}).^2+(interpolate2htoh3D(uzz{3})-uzz{4}).^2)*h^3))));

    Ep1=sqrt(sum(sum(sum(((interpolate2htoh3D(Pp{1})-Pp{2}).^2)*(4*h)^3))));
    Ep2=sqrt(sum(sum(sum(((interpolate2htoh3D(Pp{2})-Pp{3}).^2)*(2*h)^3))));
   % Ep3=sqrt(sum(sum(sum(((interpolate2htoh3D(Pp{3})-Pp{4}).^2)*h^3))));

    Ex1=sqrt(sum(sum(sum(((Xx{1}-Xx{2}).^2)*h^3))));
    Ex2=sqrt(sum(sum(sum(((Xx{2}-Xx{3}).^2)*h^3))));
    %Ex3=sqrt(sum(sum(sum(((Xx{3}-Xx{4}).^2)*h^3))));

    Ru1=log2(Eu1/Eu2);
    Ru2=log2(Eu2/Eu3);
    Ru=mean([Ru1,Ru2]);

    Rp1=log2(Ep1/Ep2);
    Rp2=log2(Ep2/Ep3);
    Rp=mean([Rp1, Rp2]);

    Rx1=log2(Ex1/Ex2);
    Rx2=log2(Ex2/Ex3);
    Rx=mean([Rx1,Rx2]);
    
    
    e0=sqrt(sum(sum(sum((uzz{1}-uz_exact{1}).^2*(8*h)^3))));
    e1=sqrt(sum(sum(sum((uzz{2}-uz_exact{2}).^2*(4*h)^3))));
    e2=sqrt(sum(sum(sum((uzz{3}-uz_exact{3}).^2*(2*h)^3))));
    e3=sqrt(sum(sum(sum((uzz{4}-uz_exact{4}).^2*(h)^3))));
    
    E0=log2(e0/e1);
    E1=log2(e1/e2);
    E2=log2(e2/e3);

end

save('ConvergenceTestSpaceBiofilm.mat')
%p=profile('info');
%save profile_results p;
%profile off
runid

