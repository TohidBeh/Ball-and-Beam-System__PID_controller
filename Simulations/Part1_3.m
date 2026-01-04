clc
clear
close all

tspan = 0 : 0.01 : 5;

u = @(t) 0.1*heaviside(t);

x0 = [0; 0; 0; 0];


[t, x] = ode45(@(t, x) myODE(x, u(t)), tspan, x0);

figure
plot(t, x(:, 1))
xlabel('time');
title("x(t)")
figure
plot(t, x(:, 2))
xlabel('time');
title("x dot(t)")
figure
plot(t, x(:, 3))
xlabel('time');
title("theta(t)")
figure
plot(t, x(:, 4))
xlabel('time');
title("theta dot(t)")

function dxdt = myODE(x, u)
    
    m = 0.27;
    r = 0.02;
    b = 1;
    k = 1e-3;
    l = 0.49;
    j_w = 14.025e-2;
    j_b = 4.32e-5;
    g = 9.81;
    
    dxdt = zeros(4, 1);
    dxdt(1) = x(2);
    dxdt(2) = ((m*r^2+j_b)*(1/r)*((2*m*x(1)*x(2)+b*l^2)*x(4) + k*l^2*x(3)-m*g*x(1)*cos(x(3))-u*l*cos(x(3)))...
        +(m*x(1)^2 +j_b+j_w)*(m*x(1)*x(4)^2+m*g*sin(x(3)))) / ((j_b/r^2 +m)*(m*x(1)^2 +j_b+j_w)- ((m*r^2 +j_b)/r)^2);
    dxdt(3) = x(4);
    dxdt(4) = ((j_b/r^2 +m)*((2*m*x(1)*x(2)+b*l^2)*x(4) + k*l^2*x(3)-m*g*x(1)*cos(x(3))-u*l*cos(x(3)))...
        +(m*r^2 +j_b)*(1/r)*(m*x(1)*x(4)^2+m*g*sin(x(3)))) / (-((j_b/r^2 +m)*(m*x(1)^2 +j_b+j_w)- ((m*r^2 +j_b)/r)^2));

end