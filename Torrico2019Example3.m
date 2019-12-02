reseta
s = tf('s');

% Example 3 ---------------------------------------------------------------
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
[A,B,C,D] = tf2ss(P1z.num{1},P1z.den{1});
A = A';
buff = B';
B = C';
C = buff;

% Pole definition
p1    = 0.95; % p1 = 0.9;
p2    = p1  ; % p2 = 0.85;
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

[G1aNum, G1aDen] = ss2tf(A,B,K,0);

nG  = G1z.num{1}; nG = nG(2:end);
nGa = G1aNum(2:end);
denV  = conv(conv([1 -beta1],[1 -beta2]),[1 -beta3]);

Ng_p1 = polyval(nG,p(1));
Ng_p2 = polyval(nG,p(2));

Nga_p1 = polyval(nGa,p(1));
Nga_p2 = polyval(nGa,p(2));

numAux = conv(conv(denV,nGa),eye(1,d+1));
denAux = nG;

numD = conv(polyder(numAux),denAux) - conv(polyder(denAux),numAux);
denD = conv(denAux,denAux);

A_v = [     1      1   1
       p(1)^3 p(1)^2 p(1)
            3      2   1];
B_v = [Kr*(1-beta1)^3
       Nga_p1*p(1)^d*(p(1)-beta1)^3/Ng_p1
       polyval(numD,p(2))/polyval(denD,p(2))];
v = A_v\B_v

V = tf([v(1) v(2) v(3) 0],denV,Ts);
% V =
%   3.123 z^3 - 6.237 z^2 + 3.114 z
%   --------------------------------
%   z^3 - 2.97 z^2 + 2.94 z - 0.9703

phiS = 0;
for i = 1:d
    phiS = phiS + K*A^(i-1)*B*z^-i;
end
phiS = minreal(phiS);
[numPhiS,denPhiS] = tfdata(phiS,'v');

numPhiS(1) = 1;
c = conv(numPhiS,denV);

A = [     1      1    1
     p(1)^3 p(1)^2  p(1)
          3      2    1];

Z = 1;
b1 = -polyval(c,Z);
Z = p(1);
b2 = -polyval(c,Z);
Z = 1;
b3 = -polyval(polyder(c),Z);

b = [b1;b2;b3];

vs = A\b;
Vs = tf([vs(1) vs(2) vs(3)],denV,Ts);
% Vs =
%   -0.1956 z^2 + 0.3752 z - 0.1796
%   --------------------------------
%   z^3 - 2.97 z^2 + 2.94 z - 0.9703

% V = Vs;

S = phiS + Vs*z^-d;
S = minreal(S);

numS = S.num{1};
denS = S.den{1};

Z = 1;
disp(polyval(numS,Z)/polyval(denS,Z))

% nS = [0 0.110050167084179 -0.434399163324191 0.643004502517625 -0.423011323349226 0.104355819571613 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00612126048742635 -0.0186592261379670 0.0189537733136556 -0.00641581016311494 0];
% dS = [1 -4.98005016708417 9.92019916332415 -9.88029850251756 4.92020017334918 -0.980050667071602 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% S  = tf(nS,dS,Ts);
% nV = [3.123305087393900 -6.237194887716494 3.113914675530637 0];
% dV = [1 -2.970000000000000 2.940300000000001 -0.970299000000001];
% V  = tf(nV,dV,Ts);
% F  = tf([0.0054 0 0],conv([1 -0.99],[1 -0.99]),Ts)/2.15;

% Call simulation
sim('Torrico2019SimuEx3')





















