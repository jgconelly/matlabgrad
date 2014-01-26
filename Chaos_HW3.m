close all; clear; format compact; clc; 

p = 0;
f = @(t,x) [x(2);-x(1)+p*(1-x(2)^2)*x(2)];
x1 = linspace(-10,10,40);
x2 = linspace(-10,10,40);
[x,y] = meshgrid(x1,x2);
u = zeros (size(x));
v = zeros (size(x));
t=0;
for i = 1:numel(x)
    xpr = f(t,[x(i);y(i)]);
    u(i)=xpr(1);
    v(i)=xpr(2);
end
quiver(x,y,u,v,'g'); figure (gcf)
xlabel ('x_1');
ylabel ('x_2');
axis tight equal

hold on
for x10 = [-2,-1,0,1,2]
for x20 = [-2,-1,0,1,2]
    [ts,xs] = ode45(f,[0,10],[x10,x20]);
    plot(xs(:,1),xs(:,2))
    plot(xs(1,1),xs(1,2),'bo') 
    plot(xs(end,1),xs(end,2),'ks')
end
end
title ('p=0; x_1(0)=-2,-1,0,1,2; x_2(0)=-2,-1,0,1,2');
hold off