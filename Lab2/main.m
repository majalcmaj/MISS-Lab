clear all
close all
graphics_toolkit ("fltk")

BOUNCES_COUNT = 5;
GRADIENT_STEP = 0.1;

G = 10;
p = [0, 0, 5];
v = [2, 2, 0];

coefs = [v(1).^2 + v(2).^2 + G / 2, 2 * (p(1) * v(1) + p(2) * v(2)) - v(3), p(1).^2 + p(2).^2 - p(3)];
delta = coefs(2).^2 - 4 * coefs(1) * coefs(3);

if delta < 0
  break;
elseif delta == 0
  impact_time = -coefs(2) / (2 * coefs(1));
else 
  impact_time = max([(-coefs(2) + sqrt(delta))/ (2 * coefs(1)), (-coefs(2) - sqrt(delta))/ (2 * coefs(1))]);
endif

impact_p = [p(1) + v(1) * impact_time, p(2) + v(2) * impact_time, p(3) + v(3) * impact_time - G * impact_time.^2 / 2];
impact_v = [v(1), v(2), v(3) - G * impact_time];

x_grad = gradient([impact_p(3), (p(1) + GRADIENT_STEP ).^2 + p(2).^2], GRADIENT_STEP);
y_grad = gradient([impact_p(3), p(1).^2 + (p(2) + GRADIENT_STEP ).^2], GRADIENT_STEP);
x_grad = [x_grad(1), 0, x_grad(2)];
y_grad = [0, y_grad(1), y_grad(2)];

norm_v = cross(x_grad, y_grad);
norm_v = norm_v ./ norm(norm_v);

[X,Y] = meshgrid(-2:0.3:2,-2:0.3:2);
Z = X.^2 + Y.^2;
hold on
s1 = surf(X,Y,Z)
set(s1,'facealpha',0.2)
norm_versor = vertcat(impact_p, impact_p + norm_v);
plot3(norm_versor(:,1), norm_versor(:,2), norm_versor(:,3), 'r', 'Linewidth', 3);
% plot3(impact_p(1), impact_p(2), impact_p(3), 'r*', 'Linewidth', 3);