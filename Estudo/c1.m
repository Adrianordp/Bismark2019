% Adriano, Aloísio, Ítalo, Judá e Nickson
reseta, format short g

s = tf('s');

% Example 3 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 60   ; % Total time
distTime   = 30   ; % Dirtubance time
distAmp    =-1  ; % Dirtubance amplitude
noiseTime  = 55   ; % Noise input time
noisePower = 1e-5 ; % Noise input power
refAmp     = 1    ; % Reference
uLim       = 90.15;

% System
P1s = (-s+1)/(2*s+1)/(3*s+1)*exp(-s);
G1s = (-s+1)/(2*s+1)/(3*s+1);

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
r1     = .975; % .975
r2     = .97;  % .97
rho    = [r1 r2];
alphaf = .95; % .95
betaf  = .90; % .90
betav  = .8; % .95
betavs = [.95 .965]; % .965
nz     = 2;
nDv    = length(betav);
nDvs   = length(betavs);
nNv    = 1;
nNvs   = 1;

% Control gain
K = acker(A,B,rho);

% Reference filter
Kr = 1/(C/(eye(size(A))+B*K-A)*B);
F  = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);

[nGa, dGa] = ss2tf(A,B,K,0);

nG  = G1z.num{1} ;
dG  = G1z.den{1} ;
dV  = [1 -betav ];
% dVs = [1 -betavs];
dVs = conv([1 -betavs(1)],[1 -betavs(2)]);

Ng_p1 = polyval(nG,p(1));
Ng_p2 = polyval(nG,p(2));

Nga_p1 = polyval(nGa,p(1));
Nga_p2 = polyval(nGa,p(2));

phi = 0;
for i = 1:d
    phi = phi + K*A^(i-1)*B*z^-i;
end
phi = minreal(phi);
% rho = [.91 .92];
x  = [rho 1 1];
Nrho = length(rho);
N  = length(x);
% N  = nv+nvs+2;
Av = zeros(N);
Bv  = zeros(N,1);
for i = 1:N
    Pis = 1
    for j = 1:nDvs
        Pis = Pis*(x(i)-betavs(j));
    end
    Pi = 1
    for j = 1:nDv
        Pi = Pi*(x(i)-betav(j));
    end
    if i <= Nrho
%          Av(i,:) = [rho(i)^(-d+1)/(rho(i)-betavs) rho(i)^-d/(rho(i)-betavs) evalfr(-Pr,rho(i))*rho(i)/(rho(i)-betav) evalfr(-Pr,rho(i))/(rho(i)-betav)];
        for j = 1:nNvs+1
            Av(i,j) = rho(i)^(nNvs-d-j+1)/Pis
        end
        for j = 1:nNv+1
            Av(i,j+nNvs+1) = -evalfr(Pr,rho(i))*rho(i)^(nNv-j+1)/Pi
        end
        Bv(i) = evalfr(phi,rho(i))+1;
    elseif i == Nrho+1
        Av(i,:) = [ones(1,nNvs+1) zeros(1,nNv+1)];
        Bv(i) = (evalfr(phi,1)+1)*Pis;
    else
        Av(i,:) = [zeros(1,nNvs+1) ones(1,nNv+1)];
        Bv(i) = Kr*Pi;
    end
end

v = Av\Bv;
% v = pinv(Av)*Bv;

Vs = tf([v(1) v(2)],dVs,Ts)
V  = tf([v(3) v(4)],dV,Ts)

S   = phi-Vs*z^-d;
S   = minreal(S);
Ceq = 1/(S+1);
YR  = series(feedback(series(Ceq,Pr),V),Kr);
zpk(YR)
if (evalfr(V,1)-Kr < 1e-10)
    disp('V(1) == Kr, ok!')
else
    disp('V(1) =/= Kr, error!')
end

if (evalfr(phi-Vs+1,1) < 1e-10)
    disp('Ceq(1) == inf, ok!')
else
    disp('Ceq(1) =/= inf, error!')
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