close all
clear all
graphics_toolkit ("fltk")
f = @(x,y) x.^2 + y.^2;
[xx,yy] = meshgrid(-5:0.25:5);
[fx,fy] = gradient(f(xx,yy),0.25);
x0 = 1;
y0 = 2;
t = (xx == x0) & (yy == y0);
indt = find(t);
fx0 = fx(indt);
fy0 = fy(indt);
z = @(x,y) f(x0,y0) + fx0*(x-x0) + fy0*(y-y0);
surf(xx,yy,f(xx,yy))
hold on
surf(xx,yy,z(xx,yy))
plot3(1,2,f(1,2),'r*')