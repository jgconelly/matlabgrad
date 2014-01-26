close all; clear; format compact;
%% Problem 1B
A = [-2 2;-2 -7];
[V,D]=eig(A)
x0 = [1;1];
[t,x] = ode45(@Prob1, [0 2], x0);
x1=x(:,1);
x2=x(:,2);
plot(t,x1,t,x2,'-');
legend('x1','x2');
title ('Problem 1b');

%% Problem 2B
clear;
A = [4 -3;6 -5];
[V,D]=eig(A)
x0 = [1;2];
[t,x] = ode45(@Prob2, [0 2], x0);
x1=x(:,1);
x2=x(:,2);
figure;
plot(t,x1,'.',t,x2,'-');
legend('x1','x2');
title ('Problem 2b');

%% Problem 2C
clear;
A = [4 -3;6 -5];
[V,D]=eig(A)
x0 = [1;2];
[t,x] = ode45(@Prob2, [0 20], x0);
x1=x(:,1);
x2=x(:,2);
figure;
plot(t,x1,'.',t,x2,'-');
legend('x1','x2');
title ('Problem 2c');

%% Problem 2D
clear;
A = [4 -3;6 -5];
[V,D]=eig(A)
x0 = [1;2];
[t,x] = ode45(@Prob2, [0 40], x0);
x1=x(:,1);
x2=x(:,2);
figure;
plot(t,x1,'.',t,x2,'-');
legend('x1','x2');
title ('Problem 2d');

%% Problem 3
clear;
A = [0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;-12 -4 15 5 -3];
[V,D]=eig(A)
x0 = [1.2500;-.50000;1;0;0];
[t,x] = ode45(@Prob3, [0 1000], x0);
x1=x(:,1);
x2=x(:,2);
x3=x(:,3);
x4=x(:,4);
x5=x(:,5);
figure;
plot(t,x1,t,x2,t,x3,t,x4,t,x5);
legend('x1','x2','x3','x4','x5');
title ('Problem 3');

%% Problem 4A
clear; 
t=0:.01:1;
x1=   exp(-2*t)   + 2*exp(-4*t) + exp(t)                ;
x2= 2*exp(-4*t)   + 2*exp(t)                            ;
x3=   exp(-2*t)   +   exp(-4*t) + exp(-0.5*t) + exp(t)  ;
x4=   exp(-2*t)   + 2*exp(t)                            ;
figure;
plot(t,x1,t,x2,t,x3,t,x4);
legend ('x1','x2','x3','x4');
title ('Problem 4a');
T = [1 2 0 1;0 2 0 0;1 1 1 1;1 0 0 0];
T1=inv(T);
lambda=[-2 0 0 0;0 -4 0 0;0 0 -.5 0;0 0 0 1];
A=T*lambda*T1

%% Problem 4B
clear;
T = [1 2 0 1;0 2 0 0;1 1 1 1;1 0 0 0];
T1=inv(T);
lambda=[-2 0 0 0;0 -4 0 0;0 0 -.5 0;0 0 0 1];
A=T*lambda*T1
[V,D] = eig(A)
D=lambda;
x0 = [0;0;1;0];
[t,x] = ode45(@Prob4, [0 2], x0);
x1=x(:,1);
x2=x(:,2);
x3=x(:,3);
x4=x(:,4);
figure;
plot(t,x1,'-',t,x2,'--',t,x3,'-',t,x4,'.');
legend('x1','x2','x3','x4');
title ('Problem 4b');





