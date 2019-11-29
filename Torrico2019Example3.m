reseta
s = tf('s');

% Example 2 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 25  ; % Total time
distTime   = 10  ; % Dis
noiseTime  = 20  ; % Noise input time
noisePower = 1e-7; % Noise input power

% System
P1s = (1-0.2*s)/s/(s-1)*exp(-0.2*s);
G1s = (1-0.2*s)/s/(s-1);

% Discretization
Ts = 0.01;
z = tf('z',Ts);
P1z = c2d(P1s,Ts);
d = P1z.InputDelay;
G1z = c2d(G1s,Ts);

% Space-state system (must be a observable cannonical form)
[A,B,C,D]  = tf2ss(P1z.num{1},P1z.den{1});
A = A';
buff = B';
B = C';
C = buff;

% Pole definition
p1    = 0.95;
p2    = p1;
beta  = [p1 p2];
alphaf= 0.00;
betaf = 0.99;
beta1 = 0.99;
beta2 = 0.99;
beta3 = 0.99;
p  = roots(P1z.den{1}); % rounded open-loop poles
nz = 2;

% Control gain
K = acker(A,B,beta); % K = [110.0255 101.6672];

% Reference filter
Kr = 1/(C/(eye(size(A))+B*K-A)*B);
F = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);

%
G1a = K*((z*eye(2)-A)\B);% g_ = K*inv(z*eye(2) - A)*B
[pza,pa,ka]=zpkdata(G1a,'v');

[G1aNum G1aDen] = ss2tf(A,B,K,0);

% % Computing V(z) (v0, v1 and v2)
% Ng_p1 = polyval(G1z.num{1},p(1));
% Ng_p2 = polyval(G1z.num{1},p(2));
% 
% Nga_p1 = polyval(G1aNum,p(1)^-1);
% Nga_p2 = polyval(G1aNum,p(2)^-1);
% 
% A_v = [1       1       1
%        1 p(1)^-1 p(1)^-2
%        1 p(2)^-1 p(2)^-2];
% B_v = [Kr*(1-beta1)^d
%        Nga_p1*p(1)^d*(1-p(1)^-1*beta1)^d/Ng_p1
%        Nga_p2*p(2)^d*(1-p(2)^-1*beta1)^d/Ng_p2];
% vs = A_v\B_v;

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
V = (v0*z^3 + v1*z^2 + v2*z)/(z-beta1)/(z-beta2)/(z-beta3)

% [r,p,k] = residue(conv(G1z.num{1},V.num{1}),conv(G1z.den{1},V.den{1}));
% VV = tf(conv(G1z.num{1},V.num{1}),conv(G1z.den{1},V.den{1}),Ts);
% NN = VV.num{1};
% DD = VV.den{1};
% Dg = G1z.den{1};
% Ng = G1z.num{1};
% Dv = conv(conv([1 -beta1],[1 -beta2]),[1 -beta3]);
% Nv = V.num{1};
% 
% b1 = polyval(NN,beta1)/polyval(Dg,beta1);
% b2 = polyval(NN,beta2)/polyval(Dg,beta2);
% b3 = polyval(NN,beta3)/polyval(Dg,beta3);
% 
% AA = [beta1^2 beta1 1
%       beta2^2 beta2 1
%       beta3^2 beta3 1];
% Nvs = (AA\[b1;b2;b3])'
% Vs  = tf(Nvs,Dv,Ts)

% S = Vs*z^-d;
% for i = 1:d
%     S = S + K*A^(i-1)*B*z^-i;
% end

% Computing S(z)
S = tf([conv(G1aNum,V.den{1}) zeros(1,d)] - [zeros(1,d) conv(V.num{1},G1z.num{1})],...
       [conv(G1z.den{1},V.den{1}) zeros(1,d)],Ts);
S = minreal(S);

% Call simulation
sim('Torrico2019SimuEx2')





















