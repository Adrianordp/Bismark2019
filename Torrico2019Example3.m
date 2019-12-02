reseta
s = tf('s');

% Example 3 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 25  ; % Total time
distTime   = 100 ; % Dis
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

[G1aNum G1aDen] = ss2tf(A,B,K,0);

nG  = G1z.num{1}; nG = nG(2:end);
nGa = G1aNum(2:end);

Ng_p1 = polyval(nG,p(1));
Ng_p2 = polyval(nG,p(2));

Nga_p1 = polyval(nGa,p(1));
Nga_p2 = polyval(nGa,p(2));

% Computing V(z) (v0, v1 and v2)
A_v = [  1      1      1
       p(1)^3 p(1)^2 p(1)
       p(2)^3 p(2)^2 p(2)];
B_v = [Kr*(1-beta1)^3
       Nga_p1*p(1)^d*(p(1)-beta1)^3/Ng_p1
       Nga_p2*p(2)^d*(p(2)-beta1)^3/Ng_p2];
vs = A_v\B_v

V = (vs(1)*z^3 + vs(2)*z^2 + vs(3)*z)/(z-beta1)/(z-beta2)/(z-beta3);
denV  = V.den{1};

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

A = [     1      1   1
     p(1)^3 p(1)^2 p(1)
          3      2   1 ];

Z = 1;
b1 = -polyval(c,Z);
Z = p(1);
b2 = -polyval(c,Z);
Z = 1;
b3 = -polyval(polyder(c),Z);

b = [b1;b2;b3];

vs = A\b;
Vs = tf([vs(1) vs(2) vs(3) 0],denV,Ts);


% Computing S(z)
S = tf([conv(G1aNum,V.den{1}) zeros(1,d)] - [zeros(1,d) conv(V.num{1},G1z.num{1})],...
       [conv(G1z.den{1},V.den{1}) zeros(1,d)],Ts);
S = minreal(S);

Vs = tf([3.123 -6.237 3.114 0],denV,Ts)
% V = Vs;
F = tf([0.0054 0 0],conv([1 -0.99],[1 -0.99]),Ts);

S = phiS + Vs*z^-d;
S = minreal(S);

numS = S.num{1};
denS = S.den{1};

Z = 1;
polyval(numS,1)/polyval(denS,1)

% Call simulation
sim('Torrico2019SimuEx3')





















