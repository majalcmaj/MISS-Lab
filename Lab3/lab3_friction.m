%Wyklad_5 strony 15-20
close all; clear all;
u = 0.01; %liniowa gestosc masy [kg / m]
T = 100; % napiecie struny [N]
b = 0.5; % Współczynnik tlumienia
L = 1; %dlugosc struny
M = 100; %liczba czesci podzialu dlugosci
N = 1000; %liczba krokow czasu
TIME = 0.1;

v = sqrt(T/u); %obliczamy predkosc struny
dx = L/M;
dt = TIME/N;
y = zeros(M,N); %rownanie struny od punktu dlugosci xm i czasu tn y(xm,tn)

% Tlumienie
beta = b / 2 / u;
uu = 1 - beta * dt;
q = 1 + beta * dt;

%uzupelniamy pierwszy rzad poczatkowymi wartosciami (struna dla t=1 <- tak na prawde zero ale Matlab)
for i = 1:M
    xi = L*(i-1)/(M-1); %indeksowanie od 1 w Matlabie
    y(i, 1) = -2*xi*xi + 0.8*xi + 0.92;
end

%uzupelniamy wartosci na brzegach struny dla kazdego t
y(1, :) = y(1,1);     %xm = 0 ??????
y(M, :) = y(M, 1); %xm = L

%teraz z wartosci poczatkowych wyliczamy strune w drugiej chwili czasu
%for i = 2:M-1
%   xi = L*(i-1)/(M-1);
%   y(i, 2) =  (v*dt/dx)^2 / 2 * ( y(i+1,1) - 2*y(i,1) + y(i-1,1)) + y(i,1) + dt * sin(2*pi*xi);
%end
p = (v * dt / dx) ^ 2;
g = zeros(M,1);
y(2:M-1, 2) = (p/2) * (y(3:M,1) - 2 * y(2:M-1, 1) + y(1:M-2, 1)) + y(2:M-1, 1) + uu * dt * g(2:M-1);

%nareszcie majac 2 pierwsze chwile czasu obliczamy y dla pozostalych t
for n = 2:N-1
  y(2:M-1,n+1) = p / q * (y(3:M, n) - 2*y(2:M-1,n) + y(1:M-2,n)) + 2 / q * y(2:M-1,n) - uu / q * y(2:M-1, n-1);
%    for m = 2:M-1
%        y(m, n+1) = 1/(2*b*dt+1)*( (dt*v/dt)^2*(y(m+1,n)-2*y(m,n)+y(m-1,n)) + 2*b*dt*y(m,n) + 2*y(m,n) - y(m,n-1) );
%    end
end

[X, T] = meshgrid(1:N, M/2:M);
size(X)
size(T)
size(y(M/2:M,:))
surf(X, T, y(M/2:M,:), 'EdgeColor','none','LineStyle','none','FaceLighting','phong')

%for t = 1:N
%  plot(y(:,t));
%  input("asd")
%end


