function parsave(w,dt,fmax,b,addlvisc,connectdist,numtimesteps,charLength)
str1=num2str(w);
str2=num2str(1/dt);
str3=num2str(fmax/1000);
str4=num2str(b/1000);
str5=num2str(addlvisc);
str6=num2str(connectdist);

if w>25
	JMain3DSimShroom2
else 
	JMain3DSimShroom2pc2
end
[G1,G2,Delta,max_stress,min_stress]=Analyze_Data2_3D(w,dt,e0,v0*visc0/charLength*vShear-eShear,fStraintop+bStraintop,numtimesteps);

runid=['w',str1,'_f',str3,'_b',str4,'_u',str5,'_cnd185_dt',str2];

save(['/lustre/janus_scratch/jast0817/',runid,'.mat']);

end
