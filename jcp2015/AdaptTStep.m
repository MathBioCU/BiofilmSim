function [speed,time,dt]= AdaptTStep(speed,time,dt,accel,c3)

speed1=speed(c3-1)+dt(c3)*accel;
time1=time(c3-1)+dt(c3);


t4=time(c3-4); s4=speed(c3-4);
t3=time(c3-3); s3=speed(c3-3);
t2=time(c3-2); s2=speed(c3-2);
t1=time(c3-1); s1=speed(c3-1);

L1=@(t) (t-t2)/(t1-t2)*(t-t3)/(t1-t3)*(t-t4)/(t1-t4);
L2=@(t) (t-t1)/(t2-t1)*(t-t3)/(t2-t3)*(t-t4)/(t2-t4);
L3=@(t) (t-t1)/(t3-t1)*(t-t2)/(t3-t2)*(t-t4)/(t3-t4);
L4=@(t) (t-t1)/(t4-t1)*(t-t2)/(t4-t2)*(t-t3)/(t4-t3);

p=@(t) L1(t)*s1+L2(t)*s2+L3(t)*s3+L4(t)*s4;

tol=0.005*speed(c3-1);
dt(c3)=dt(c3)*0.9*tol/(abs(p(time1)-speed1));
time(c3)=time(c3-1)+dt(c3);

speed(c3)=speed(c3-1)+dt(c3)*accel;
end

