close all; clear all;
mu=0.01; % liniowa gestosc masy
tau=100; % napiecie struny
L=1; % dlugosc struny
M=100; % odcinki struny
N=2000; % odcinki czasu
x=linspace(0,L,M)';
v=sqrt(tau/mu);
dx=x(2)-x(1);
dt=dx/v;
p=(v*dt/dx)^2;
s=-2 * x.^2 + 0.8 * x + 0.92; 
g= sin(2*pi*x); % tutaj warunek co do pochodnej z y po czasie dla t=0
h= zeros(1,N) + 0.8; % tutaj dodac kat(jesli jest staly)
r= zeros(1,N) - 0.28; % tutaj dodac polozenie przy x=0 (stale)

% Tlumienie
b = 0.5;
beta = b/2/mu;
q = 1 + beta * dt;
u = 1 - beta * dt;

f(:,1)=s;
f(1,2)=p*(f(2,1)- f(1,1)-dx*h(1)) + f(1,1) +dt*g(1);
f(2:M-1,2)=p/2*(f(3:M,1)-2*f(2:M-1,1)+ f(1:M-2,1))+f(2:M-1,1)+dt*g(2:M-1);
f(M,2)=r(2);
for n=2:N-1
    f(1,n+1)=2*p*(f(2,n)-f(1,n)- dx*h(n) )+2*f(1,n)-f(1,n-1);
    f(2:M-1,n+1)=p/q*(f(3:M,n)-2*f(2:M-1,n)+f(1:M-2,n))+2/q*f(2:M-1,n)-u/q*f(2:M-1,n-1);
    f(M,n+1)=r(n+1);
end

[X, T] = meshgrid(1:N, 1:M);
size(X)
size(T)
size(f(1:M,:))
surf(X, T, f(1:M,:), 'EdgeColor','none','LineStyle','none','FaceLighting','phong')

eq_x = [0,1];
eq_y_analytic = 0.8 * eq_x - 0.28
eq_y_computed = [ f(M,N) f(1,N) ]
plot(eq_y_analytic);
hold on;
plot(eq_y_computed);