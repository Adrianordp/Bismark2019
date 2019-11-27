reseta
s = tf('s');

% Example 2 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 30     ; % Total time
noiseTime  = 99     ; % Noise input time
noisePower = 0.00001; % Noise input power

% System
P1s = 2/(3*s-1)/(s-1)*exp(-0.3*s);
G1s = 2/(3*s-1)/(s-1);

% Discretization
Ts = 0.05;
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
p1    = exp((-1/0.85)*Ts);
p2    = p1;
beta  = [p1 p2];
alphaf= 0.94;
betaf = 0.8985;
beta1 = 0.529;
beta2 = 0.53;
beta3 = 0.531;
p  = round(roots(P1z.den{1})*1e8)/1e8; % rounded open-loop poles
% p  = roots(P1z.den{1}); % rounded open-loop poles
nz = 2;

% Control gain
K = acker(A,B,beta); % K = [110.0255 101.6672];

% Reference filter
Kr = 1/(C/(eye(size(A))+B*K-A)*B);
F = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);

%
G1a = K*((z*eye(2)-A)\B);% g_ = K*inv(z*eye(2) - A)*B
[pza,pa,ka]=zpkdata(G1a,'v');

G1aNum = conv(ka,[1 -pza(3)]);

Ng_p1 = polyval(G1z.num{1},p(1));
Ng_p2 = polyval(G1z.num{1},p(2));

% Nga_p1 = polyval(G1a.num{1},p(1));
% Nga_p2 = polyval(G1a.num{1},p(2));

Nga_p1 = polyval(G1aNum,p(1)^-1);
Nga_p2 = polyval(G1aNum,p(2)^-1);

A_v = [1       1       1
       1 p(1)^-1 p(1)^-2
       1 p(2)^-1 p(2)^-2];
B_v = [Kr*(1-beta1)^d
       Nga_p1*p(1)^d*(1-p(1)^-1*beta1)^d/Ng_p1
       Nga_p2*p(2)^d*(1-p(2)^-1*beta1)^d/Ng_p2];
vs = A_v\B_v;
%

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
% V = (vs(1)*z^3 + vs(2)*z^2 + vs(3)*z)/(z-beta1)/(z-beta2)/(z-beta3);
% V = (231.6095*z^3 - 455.3408*z^2 + 223.9279*z)/(z-beta1)^3;

[r,p,k] = residue(conv(G1z.num{1},V.num{1}),conv(G1z.den{1},V.den{1}));
VV = tf(conv(G1z.num{1},V.num{1}),conv(G1z.den{1},V.den{1}),Ts);
NN = VV.num{1};
DD = VV.den{1};
Dg = G1z.den{1};
Ng = G1z.num{1};
Dv = conv(conv([1 -beta1],[1 -beta2]),[1 -beta3]);
Nv = V.num{1};

b1 = polyval(NN,beta1)/polyval(Dg,beta1);
b2 = polyval(NN,beta2)/polyval(Dg,beta2);
b3 = polyval(NN,beta3)/polyval(Dg,beta3);

% AA = [1 beta1^-1 beta1^-2
%       1 beta2^-1 beta2^-2
%       1 beta3^-1 beta3^-2];
AA = [beta1^2 beta1 1
      beta2^2 beta2 1
      beta3^2 beta3 1];
Nvs = (AA\[b1;b2;b3])'
Vs = tf(Nvs,Dv,Ts)

% eV = eig(V);
% pos = [3;4;5];
% 
% Vs = 0;
% for i = 1:length(pos)
%     Vs = Vs + r(pos(i))/(z-p(pos(i)));
% end

% V = Vs;
% Vs = -(0.080908)*(z-0.1255)*(z-2.467)/(z-0.53)^3;
% Vs = -(0.087574)*(z-4.108)*(z-0.1098)/(z-0.53)^3;
% V = Vs;
S = Vs*z^-d;
for i = 1:d
    S = S + K*A^(i-1)*B*z^-i;
end
S = minreal(S);
% S = (conv(G1aNum,Vs.den{1}) - conv(Vs.num{1},G1z.num{1})*z^-d)/conv(G1z.den{1},Vs.den{1})
ns=[0 0.182331139052535 -0.0927609385675930 0.0530004141616560 0.0299933867604861 0.0321617711957800 0.0344519289926015 -0.160485363778323 -0.182515337817140]
ds=[1 -1.59 0.8427 -0.148877 0 0 0 0 0]

% ng_ = [0    0.1823   -0.1799];
% dV  = [1.0000   -1.5900    0.8427   -0.1489]
% dn  = 6;
% dG  = [1.0000   -2.0681    1.0689]

Snum = [conv(ng_,dV) zeros(1,dn)] - [zeros(1,dn) conv(nV,nG)]
Sden = [conv(dG,dV) zeros(1,dn)];
S = minreal(tf(Snum,Sden,Ts))

S = tf(ns,ds,Ts);
sim('Torrico2019SimuEx2')





















