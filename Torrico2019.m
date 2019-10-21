reseta

Tsim = 60;

Ts = 0.1;
s = tf('s');
z = tf('z',Ts);
P1s = (-s+1)/(2*s+1)/(3*s+1)*exp(-s);
step(P1s,Tsim)

P1z = c2d(P1s,Ts);
d = P1z.InputDelay;

[A,B,C,D]  = tf2ss(P1z.num{1},P1z.den{1});
A = A';
buff = B';
B = C';
C = buff;

beta = [0.92 0.92];

betaf = 0.89;
beta1 = 0.9;
beta2 = 0.9;
beta3 = 0.74;
p1 = 0.92;
p2 = 0.92;

K = acker(A,B,beta)
Kr = inv(C/(eye(size(A))+B*K-A)*B);
F = Kr*(1-betaf)^2*z^2/(z-betaf)^2;


% -------------

V = (5.155*z^3 - 9.848*z^2 + 4.704*z)/(z-beta1)^2/(z-beta3);

sim('Torrico2019Simu')

syms v0 v1 v2
vv = [v0 v1 v2]';
V1 = (v0 + v1 + v2)/(1-beta1)/(1-beta2)/(1-beta3);
Vp1 = (v0 + v1*p1^-1 + v2*p1^-2)/(1-beta1)/(1-beta2)/(1-beta3);
Vp2 = (v0 + v1*p2^-1 + v2*p2^-2)/(1-beta1)/(1-beta2)/(1-beta3);

S1 = 1 + (K - V1*C)*((eye(size(A)) - A)\B);
Sp1 = 1 + (K - Vp1*C)*((eye(size(A)) - A)\B);
Sp2 = 1 + (K - Vp2*C)*((eye(size(A)) - A)\B);