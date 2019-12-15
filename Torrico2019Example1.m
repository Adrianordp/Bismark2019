reseta
s = tf('s');

% Example 1 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 60  ; % Total time
noiseTime  = 55  ; % Noise input time
noisePower = 1e-5; % Noise input power

% System
P1s = (-s+1)/(2*s+1)/(3*s+1)*exp(-s);

% Discretization
Ts  = 0.1;
z   = tf('z',Ts);
P1z = c2d(P1s,Ts);
d   = P1z.InputDelay;

% Space-state system (must be a observable cannonical form)
[A,B,C,D]  = tf2ss(P1z.num{1},P1z.den{1});
A = A';
buff = B';
B = C';
C = buff;

% Pole definition
beta   = [0.91 0.92]; % Desirable poles for closed loop system
alphaf = 0.00       ; % Reference filter (zero)
betaf  = 0.89       ; % Reference filter (pole)
beta1  = 0.90       ; % Robustness filter tunning
beta2  = 0.90       ; % Robustness filter tunning
beta3  = 0.74       ; % Robustness filter tunning
p      = round(roots(P1z.den{1})*1e8)/1e8; % rounded open-loop poles
% p      = roots(P1z.den{1}); % rounded open-loop poles
nz     = 2          ; % Order of the alphaf filter for F(z)

% Control gain
K = acker(A,B,beta);

% Reference filter
Kr = inv(C/(eye(size(A))+B*K-A)*B);
F  = Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz;

% Solution for V(z)
syms v0 v1 v2
vv  = [v0 v1 v2]';
% 1+S = 0, for z = 1
Z  = 1;
V1 = (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);
S1 = 1 + (K - Z^-d*V1 *C)*((eye(size(A))*Z - A)\B);
% 1+S = 0, for z = pole1
Z  = p(1);
Vp1= (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);
Sp1= 1 + (K - Z^-d*Vp1*C)*((eye(size(A))*Z - A)\B);
% 1+S = 0, for z = pole2
Z  = p(2);
Vp2= (v0 + v1*Z^-1 + v2*Z^-2)/(1-beta1*Z^-1)/(1-beta2*Z^-1)/(1-beta3*Z^-1);
Sp2= 1 + (K - Z^-d*Vp2*C)*((eye(size(A))*Z - A)\B);
% Solve linear problem
V  = solve(S1==0,Sp1==0,Sp2==0,v0,v1,v2);
v0 = eval(V.v0);
v1 = eval(V.v1);
v2 = eval(V.v2);
V  = (v0*z^3 + v1*z^2 + v2*z)/(z-beta1)/(z-beta2)/(z-beta3);

% Simulate curves
sim('Torrico2019Simu')