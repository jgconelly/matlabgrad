close all; clear; 
t = -10:.01:10;
p = [-10:.01:10];
J1 = -2*t;
J2 = -p;
z = 0;
plot(t,J1,'g',t,J2,'r',t,z,'-.b')
legend('-2*t','-p','zero')
