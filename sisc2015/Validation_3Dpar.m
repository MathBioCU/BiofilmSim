clear all 
fr=[0.1254 0.4991 1.578 4.991 7.91 9.958 12.54 19.87 25.01 31.49 39.64 49.91 62.83];
dt=1/800;
numtimesteps=800;
charLength=10*10^-6;
fmax=380000;
addlvisc=500;
b=20;
connectdist=3/18.5;
parpool

parfor i=3:13
    parsave(fr(i),dt,fmax,b,addlvisc,connectdist,numtimesteps,charLength);
end
Complete='Complete'
delete(gcp);

