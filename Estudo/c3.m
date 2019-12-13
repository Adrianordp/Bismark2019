% Adriano, Aloísio, Ítalo, Judá e Nickson
reseta, format short g

s = tf('s');

% Example 3 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 100 ; % Total time
distTime   = 40 ; % Dirtubance time
distAmp    =-0.2; % Dirtubance amplitude
noiseTime  = 80  ; % Noise input time
noisePower = 1e-5; % Noise input power
refAmp     = 1   ; % Reference
uLim       = 90.15;

% System
% P1s = (-s+1)/(2*s+1)/(3*s+1)*exp(-s);
% G1s = (-s+1)/(2*s+1)/(3*s+1);
P1s = (1)/(2*s+1)/(3*s+1)*exp(-s);
G1s = (1)/(2*s+1)/(3*s+1);
% P1s = (1-0.2*s)/s/(s-1)*exp(-0.2*s);
% G1s = (1-0.2*s)/s/(s-1);

% Discretization
Ts  = 0.1          ;
z   = tf('z',Ts)   ;
Pr  = c2d(P1s,Ts)  ;
d   = Pr.InputDelay;
G1z = c2d(G1s,Ts)  ;

% Space-state system (must be a observable cannonical form)
[A,B,C,D] = tf2ss(Pr.num{1},Pr.den{1});
A    = A'  ;
buff = B'  ;
B    = C'  ;
C    = buff;

Gss = ss(A,B,C,D,Ts);

% Pole definition
p      = roots(Pr.den{1}); % open-loop poles
ze     = roots(Pr.num{1});
r1     = 0.975;
r2     = 0.970;
rho    = [r1 r2];
alphaf = 0.000;
betaf  = 0.9;
betav  = .9;
betasv = .95;
nz     = 2;

% Control gain
K = acker(A,B,rho); % K = [110.0255 101.6672];
% Q = eye(2); R = .1;
% K = lqr(Gss,Q,R)

% Reference filter
Kr = 1/(C/(eye(size(A))+B*K-A)*B);
F  = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);

[nGa, dGa] = ss2tf(A,B,K,0);

nG  = G1z.num{1} ;
dG  = G1z.den{1} ;
dV  = [1 -betav ];
dVs = [1 -betasv];

Ng_p1 = polyval(nG,p(1));
Ng_p2 = polyval(nG,p(2));

Nga_p1 = polyval(nGa,p(1));
Nga_p2 = polyval(nGa,p(2));

phi = 0;
for i = 1:d
    phi = phi + K*A^(i-1)*B*z^-i;
end
phi = minreal(phi);
% rho = [.985 .965];
x  = [rho 1 1];
Nrho = length(rho);
N  = length(x);
Av = zeros(N);
Bv  = zeros(N,1);
for i = 1:N
    if i <= Nrho
        Av(i,:) = [rho(i)^(-d+1)/(rho(i)-betasv) rho(i)^-d/(rho(i)-betasv) evalfr(-Pr,rho(i))*rho(i)/(rho(i)-betav) evalfr(-Pr,rho(i))/(rho(i)-betav)];
        Bv(i) = evalfr(phi,rho(i))+1;
    elseif i == Nrho+1
        Av(i,:) = [1 1 0 0];
        Bv(i) = (evalfr(phi,1)+1)*(1-betasv);
    else
        Av(i,:) = [0 0 1 1];
        Bv(i) = Kr*(1-betav);
    end
end

v = Av\Bv;

Vs = tf([v(1) v(2)],dVs,Ts)
V  = tf([v(3) v(4)],dV,Ts)

S   = phi-Vs*z^-d;
S = minreal(S);
Ceq = V/(S+1);
YR  = minreal(Pr*Kr/(Pr*V+S+1));
zpk(YR)
if (evalfr(V,1)-Kr < 1e-10)
    disp('V(1) == Kr, ok!')
else
    disp('V(1) =/= Kr, error!')
end

if (evalfr(phi-Vs+1,1) < 1e-10)
    disp('V(1) == Kr, ok!')
else
    disp('V(1) =/= Kr, error!')
end

[zz,pp,kk] = zpkdata(YR,'v');
for i = 1:Nrho
    if (min(pp-rho(i)) < 1e-4)
        disp(['pole ' num2str(i) ' ok!'])
    else
        disp(['pole ' num2str(i) ' error!'])
    end
end

disp(max(abs(pp)))
if(max(abs(pp))>=1)
    disp('There is a pole outside the unitary circle!')
else
    % Call simulation
    sim('simuC3')
end