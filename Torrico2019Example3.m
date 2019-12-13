reseta
s = tf('s');

% Example 3 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 40  ; % Total time
distTime   = 20  ; % Dirtubance time
distAmp    =-0.06; % Dirtubance amplitude
noiseTime  = 35  ; % Noise input time
noisePower = 1e-8; % Noise input power
refAmp     = 0.2 ; % Reference
uLim       = 90.15;

% System
P1s = (1-0.2*s)/s/(s-1)*exp(-0.2*s);
G1s = (1-0.2*s)/s/(s-1);

% Discretization
Ts  = 0.01          ;
z   = tf('z',Ts)    ;
P1z = c2d(P1s,Ts)   ;
d   = P1z.InputDelay;
G1z = c2d(G1s,Ts)   ;

% Space-state system (must be a observable cannonical form)
[A,B,C,D] = tf2ss(P1z.num{1},P1z.den{1});
A    = A'  ;
buff = B'  ;
B    = C'  ;
C    = buff;

Gss = ss(A,B,C,D,Ts);

% Pole definition
p1    = 0.98; % p1 = 0.9;
p2    = p1  ; % p2 = 0.85;
beta  = [p1 p2];
alphaf= 0.00;
betaf = 0.994; % 0.99 no artigo, 0.94 real
beta1 = 0.970;
beta2 = 0.970;
beta3 = 0.970;
p     = roots(P1z.den{1}); % rounded open-loop poles
nz    = 2;

% Control gain
K = acker(A,B,beta); % K = [110.0255 101.6672];
% Q = eye(2); R = .1;
% K = lqr(Gss,Q,R)

% Reference filter
Kr = 1/(C/(eye(size(A))+B*K-A)*B);
F = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);

[nGa, dGa] = ss2tf(A,B,K,0);

nG    = G1z.num{1}   ;
dG    = G1z.den{1}   ;
dV  = conv(conv([1 -beta1],[1 -beta2]),[1 -beta3]);

Ng_p1 = polyval(nG,p(1));
Ng_p2 = polyval(nG,p(2));

Nga_p1 = polyval(nGa,p(1));
Nga_p2 = polyval(nGa,p(2));

nAux = conv(conv(dV,nGa+dG),eye(1,d+1));
dAux = nG;

numD = conv(polyder(nAux),dAux) - conv([0 polyder(dAux)],nAux);
denD = conv(dAux,dAux);

A_v = [       1      1    1
         p(1)^3 p(1)^2  p(1)
       3*p(2)^2  2*p(2)   1];
B_v = [Kr*(1-beta1)*(1-beta2)*(1-beta3)
       Nga_p1*p(1)^d*(p(1)-beta1)*(p(1)-beta2)*(p(1)-beta3)/Ng_p1
       polyval(numD,p(2))/polyval(denD,p(2))];
v = A_v\B_v;

V = tf([v(1) v(2) v(3) 0],dV,Ts);
nV = V.num{1};

phiS = 0;
for i = 1:d
    phiS = phiS + K*A^(i-1)*B*z^-i;
end
phiS = minreal(phiS);
% x = [1 p1 beta1 rand(1,2)];
x = rand(1,5);
Aa = zeros(5);
b  = zeros(5,1);
for i = 1:5
    dGZ = polyval(dG,x(i));
    dVZ = polyval(dV,x(i));
    nGZ = polyval(nG,x(i));
    nVZ = polyval(nV,x(i));
    Aa(i,:) = [[x(i)^2 x(i) 1]*dGZ [x(i) 1]*dVZ];
    b(i) = nGZ*nVZ;
end

format long g
vs = Aa\b

Vs = tf([vs(1) vs(2) vs(3)],dV,Ts)

S = phiS - Vs*z^-d; % ???
S = minreal(S);

numS = S.num{1};
denS = S.den{1};

[zz,pp,kk] = zpkdata(S,'v');
[nSt,dSt] = zp2tf(zz(1:21),pp,kk);
% [nSt,dSt] = zp2tf(zz(end-20:end),pp,kk);
St = tf(nSt,dSt,Ts);
Sb = S;
S = St;

% Call simulation
sim('Torrico2019SimuEx3')





















