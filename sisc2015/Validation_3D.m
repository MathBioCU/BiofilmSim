% clc
clear all
fr=[0.1254 0.4991 1.254 4.991 12.54 25.01 39.64 49.91 62.83];
for i=6:6%length(fr)
    
    w=fr(i);
    dt=1/100;
    numtimesteps=1;
    fmax=380000;
    addlvisc=500;
    b=20;
    connectdist=3/18.5;
    JMain3DSimShroom2
    [G1,G2,Delta,max_stress,min_stress]=Analyze_Data2_3D(w,dt,e0,v0*visc0/charLength*vShear-eShear,fStraintop+bStraintop,numtimesteps);
     str1=num2str(i);
     runid=['w',str1,'_f380_b20_u500_cd185_dt1000'];
     save(['/home/jstotsky/scratch/',runid,'.mat'])
end


