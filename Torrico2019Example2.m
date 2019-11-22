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
G1s = 2/(3*s-1)/(s-1);
G1z = c2d(G1s,Ts);

[A,B,C,D]  = tf2ss(P1z.num{1},P1z.den{1});
A = A';
buff = B';
B = C';
C = buff;

p1 = exp((-1/0.85)*Ts);
p2 = p1;
beta = [p1 p2];

alphaf= 0.94;
betaf = 0.8985;
beta1 = 0.53;
beta2 = 0.53;
beta3 = 0.53;
% p  = round(roots(P1z.den{1})*1e8)/1e8; % rounded open-loop poles
p  = roots(P1z.den{1}); % rounded open-loop poles
nz = 2;

K = acker(A,B,beta)
% K = [110.0255 101.6672];
Kr = inv(C/(eye(size(A))+B*K-A)*B);
F = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);

Ng_p1 = G1z.num{1}(1)*p1+G1z.num{1}(2)*p1+G1z.num{1}(3)*p1;
Ng_p2 = G1z.num{1}(1)*p2+G1z.num{1}(2)*p2+G1z.num{1}(3)*p2;

% Falta o Nga

A_v = [1      1       1
       1 p1^(-1) p1^(-2)
       1 p2^(-1) p2^(-2)];
B_v = [Kr*(1+beta1)^d
       Nga_p1*p1^d*(1-p1^(-1)*beta1)^d/Ng_p1
       Nga_p2*p2^d*(1-p2^(-1)*beta1)^d/Ng_p2];
vs = A_v\B_v






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
V  = solve(S1==0,Sp1==0,Sp2==0,v0,v1,v2);
v0 = eval(V.v0);
v1 = eval(V.v1);
v2 = eval(V.v2);
V = (v0*z^3 + v1*z^2 + v2*z)/(z-beta1)/(z-beta2)/(z-beta3);
% V = (231.6095*z^3 - 455.3408*z^2 + 223.9279*z)/(z-beta1)^3;

[r,p,k] = residue(conv(G1z.num{1},V.num{1}),conv(G1z.den{1},V.den{1}));

eV = eig(V);
pos = [3;4;5];

Vs = 0;
for i = 1:length(pos)
    Vs = Vs + r(pos(i))/(z-p(pos(i)))^i;
end
% V = Vs;
% Vs = -(0.080908)*(z-0.1255)*(z-2.467)/(z-0.53)^3;
Vs = -(0.087574)*(z-4.108)*(z-0.1098)/(z-0.53)^3;
% V = Vs;
S = Vs*z^-d;
for i = 1:d
    S = S + K*A^(i-1)*B*z^-i;
end
sim('Torrico2019SimuEx2')





















