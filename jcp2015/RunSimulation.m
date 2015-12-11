clc
% clear all 
fr=49.91;
%profile on
clear;
clc;
for i=1:1
    w=39.64;
    dt=1/550;
    dt=dt/w;
    dx=1/32;
    numtimesteps=5000;
    charLength=10*10^-6;
    fmax=28000;
    addlvisc=150;
    b=0.2;  %b should be fairly small before there is a significant difference i.e. 0.001. 
    B01=b;
    E=0;
    E01=E;
    sigma0=0.0;
    connectdist=0.162;
    JMain3DSim7
    [G1,G2,Delta,max_stress,min_stress]=Analyze_Data2_3D(w,w*dt,e0,v0*visc0/charLength*vShear-eShear,fStrain,numtimesteps);
    str1=num2str(w);
    str2=num2str(fmax/1000);
    str3=num2str(B01);
    str4=num2str(1/(w*dt));
    str5=num2str(connectdist);
    str6=num2str(addlvisc);
    str7=num2str(sigma0);
    str8=datestr(datetime,'mm-dd-yy');
    runid=['w',str1,'_f',str2,'_b',str3,'_cnd',str5,'_visc',str6,'_dt',str4,'_sigma',str7,'wSAR2_test',str8];
    save([runid,'.mat'])
    Complete='Complete';
    G1s=num2str(G1);
    G2s=num2str(G2);
    maxstrain=num2str(max(fStraintop+bStraintop));
    fprintf('G1 %s,G2  %s,Max Strain %s, %s\n',G1s,G2s,maxstrain,Complete);	

end
%p=profile('info');
%save profile_results p;
%profile off
runid

