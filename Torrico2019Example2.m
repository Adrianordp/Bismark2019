reseta
s = tf('s');

% Example 2 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 30     ; % Total time
noiseTime  = 99     ; % Noise input time
noisePower = 0.00001; % Noise input power

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
p  = round(roots(P1z.den{1})*1e8)/1e8; % rounded open-loop poles
nz = 2;

K = acker(A,B,beta)
Kr = inv(C/(eye(size(A))+B*K-A)*B);
F = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);

% Solution for V(z)
syms v0 v1 v2
vv  = [v0 v1 v2]';
% 1+S = 0, for z = 1
Z   = 1;
V1  = (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);
S1  = 1 + (K - Z^-d*V1 *C)*((eye(size(A))*Z - A)\B);
% 1+S = 0, for z = pole1
Z   = p(1);
Vp1 = (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);
Sp1 = 1 + (K - Z^-d*Vp1*C)*((eye(size(A))*Z - A)\B);
% 1+S = 0, for z = pole2
Z   = p(2);
Vp2 = (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);
Sp2 = 1 + (K - Z^-d*Vp2*C)*((eye(size(A))*Z - A)\B);
% Solve linear problem
V = solve(S1==0,Sp1==0,Sp2==0,v0,v1,v2);
v0 = eval(V.v0);
v1 = eval(V.v1);
v2 = eval(V.v2);
V = (v0*z^3 + v1*z^2 + v2*z)/(z-beta1)/(z-beta2)/(z-beta3);
V
% V = (231.6095*z^3 - 455.3408*z^2 + 223.9279*z)/(z-beta1)^3;
sim('Torrico2019Simu')


