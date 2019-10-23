reseta
s = tf('s');

% Example 2 ---------------------------------------------------------------
Tsim = 25;
Ts = 0.05;
z = tf('z',Ts);
P1s = 2*exp(-0.3*s)/(3*s-1)/(s-1);
P1z = c2d(P1s,Ts)
d = P1z.InputDelay;

[A,B,C,D]  = tf2ss(P1z.num{1},P1z.den{1});
A = A';
buff = B';
B = C';
C = buff;

beta = [0.942873148162653 0.942873148162653];

alphaf= 0.94;
betaf = 0.8985;
beta1 = 0.53;
beta2 = 0.53;
beta3 = 0.53;
delta = 0.0000000001;
id = 1+delta;
p1 = 0.92;
p2 = 0.92;
nz = 2;

K = acker(A,B,beta)
Kr = inv(C/(eye(size(A))+B*K-A)*B);
F = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);
V = (231.6095*z^3 - 455.3408*z^2 + 223.9279*z)/(z-beta1)^3;
% sim('Torrico2019Simu')

syms v0 v1 v2
vv = [v0 v1 v2]';
Z = 1;
V1 = (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);
Z = p1;
Vp1 = (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);
Z = id;
Vd1 = (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);

% V1 = (v0 + v1 + v2)/(1-beta1)/(1-beta2)/(1-beta3);
% Vp1 = (v0 + v1*p1^-1 + v2*p1^-2)/(1-beta1)/(1-beta2)/(1-beta3);
% Vd1 = (v0 + v1*id^-1 + v2*id^-2)/(1-beta1)/(1-beta2)/(1-beta3);

S1 = 1 + (K - V1*C)*((eye(size(A)) - A)\B)
Sp1 = 1 + (K - p1^-d*Vp1*C)*((eye(size(A))*p1 - A)\B);
Sd1 =  1 + (K - id^-d*Vd1*C)*((eye(size(A))*id - A)\B);
Sd  = (Sd1-S1)/delta

V = solve(S1,Sp1,Sd,v0,v1,v2);
v0 = eval(V.v0);
v1 = eval(V.v1);
v2 = eval(V.v2);
V = (v0*z^3 + v1*z^2 + v2*z)/(z-beta1)^3
sim('Torrico2019Simu')


