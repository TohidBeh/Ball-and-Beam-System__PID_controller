clc
clear
close all

syms x1 x2 x3 x4 t u
m = 0.27;
r = 0.02;
b = 1;
k = 1e-3;
l = 0.49;
j_w = 14.025e-2;
j_b = 4.32e-5;
g = 9.81;
threshold = 1e-4;

f1 = x2;
f2 = ((m*r^2+j_b)*(1/r)*((2*m*x1*x2+b*l^2)*x4 + k*l^2*x3-m*g*x1*cos(x3)-u*l*cos(x3))...
        +(m*x1^2 +j_b+j_w)*(m*x1*x4^2+m*g*sin(x3))) / ((j_b/r^2 +m)*(m*x1^2 +j_b+j_w)- ((m*r^2 +j_b)/r)^2);
f3 = x4;
f4 = ((j_b/r^2 +m)*((2*m*x1*x2+b*l^2)*x4 + k*l^2*x3-m*g*x1*cos(x3)-u*l*cos(x3))...
        +(m*r^2 +j_b)*(1/r)*(m*x1*x4^2+m*g*sin(x3))) / (-((j_b/r^2 +m)*(m*x1^2 +j_b+j_w)- ((m*r^2 +j_b)/r)^2));

x1_0 = 0;
x2_0 = 0;
x3_0 = 0;
x4_0 = 0;
u_0 = 0;

A = jacobian([f1; f2; f3; f4], [x1; x2; x3; x4]);
A = subs(A, [x1; x2; x3; x4; u], [x1_0; x2_0; x3_0; x4_0; u_0]);

B = jacobian([f1; f2; f3; f4], u);
B = subs(B, [x1; x2; x3; x4; u], [x1_0; x2_0; x3_0; x4_0; u_0]);


A = double(A);
B = double(B);

disp('A:');
disp(A);
disp('B:');
disp(B);


C= [1, 0, 0, 0 ; 0, 0, 1, 0]
D = [0;0]


[b, a] = ss2tf(A, B, C, D);

b1 = b(1,:);
b1 = b1.*(abs(b1)>threshold);
b2 = b(2,:);
b2 = b2.*(abs(b2)>threshold);
a = a.*(abs(a)> threshold);
x_s = tf(b1, a)
theta_s = tf(b2, a)

G = tf([-0.06993 0 24.5], [1 1.713+6.39*3.496  0.7577+58.78*3.496  0 -132.4]);

% k > 0
figure;
rlocus(G)
title("K > 0")

% k < 0
figure;
rlocus(-G)
title("K < 0")
