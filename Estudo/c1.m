% Adriano, Aloísio, Ítalo, Judá e Nickson
reseta, format short g

s = tf('s');

% Example 3 ---------------------------------------------------------------
% Simulation parameters
Tsim       = 60   ; % Total time
distTime   = 30   ; % Dirtubance time
distAmp    =-1    ; % Dirtubance amplitude
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
p      = roots(Pr.den{1})'; % open-loop poles
ze     = roots(Pr.num{1});
r1     = .975; % .975
r2     = .97; % .97
rho    = [r1 r2];
alphaf = .95; % .95
betaf  = .90; % .90
betav  = [.964 .965 .966]; % .95
betavs = [.964 .965 .966]; % .965
% betavs = [.9 .9 .74]; % .965
Nrho   = length(rho);
Np     = length(p);
nz     = 2;
nDv    = length(betav);
nDvs   = length(betavs);
nNv    = 3;
nNvs   = 2;

% Control gain
K = acker(A,B,rho);

% Reference filter
Kr = 1/(C/(eye(size(A))+B*K-A)*B);
F  = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);

[nGa, dGa] = ss2tf(A,B,K,0);

nG  = G1z.num{1};
dG  = G1z.den{1};

denV = 1;
for i = 1:nDv
    denV = conv(denV, [1 -betav(i)]);
end

denVs = 1;
for i = 1:nDvs
    denVs = conv(denVs, [1 -betavs(i)]);
end

Ng_p1 = polyval(nG,p(1));
Ng_p2 = polyval(nG,p(2));

Nga_p1 = polyval(nGa,p(1));
Nga_p2 = polyval(nGa,p(2));

phi = 0;
for i = 1:d
    phi = phi + K*A^(i-1)*B*z^-i;
end
phi = minreal(phi);
Nr  = nNvs+nNv-Nrho-Np;
if Nr < 0
    error('Aumentar ordem dos polonomios V(z) e V*(z).')
end
x  = [rho p 1 1 rand(1,Nr)/2+0.1];
x(end-1:end) = betavs(1:2);
N1 = Nrho+Np;
N2 = N1+2;
Nx = length(x);
Av = zeros(Nx,nNv+nNvs+2);
Bv = zeros(Nx,1);

Pis = 1;
for i = 1:nDvs
    Pis = conv(Pis,[1 -betavs(i)]);
end
Pis = tf(Pis,1,Ts);

Pi = 1;
for i = 1:nDv
    Pi = conv(Pi,[1 -betav(i)]);
end
Pi = tf(Pi,1,Ts);

for i = 1:Nx
    if i <= Nrho
        for j = 1:nNvs+1
            Av(i,j) = x(i)^(nNvs+1-d-j)/evalfr(Pis,x(i));
        end
        for j = 1:nNv+1
            Av(i,j+nNvs+1) = -evalfr(Pr,x(i))*x(i)^(nNv-j+1)/evalfr(Pi,x(i));
        end
        Bv(i) = evalfr(phi,x(i))+1;
    elseif i <= N1
        for j = 1:nNvs+1
            Av(i,j) = x(i)^(nNvs+1-d-j);
        end
        for j = 1:nNv+1
            Av(i,j+nNvs+1) = 0;
        end
        Bv(i) = (evalfr(phi,x(i))+1)*evalfr(Pis,x(i));
    elseif i == N1+1
        Av(i,:) = [ones(1,nNvs+1) zeros(1,nNv+1)];
        Bv(i) = (evalfr(phi,1)+1)*evalfr(Pis,x(i));
    elseif i == N1+2
        Av(i,:) = [zeros(1,nNvs+1) ones(1,nNv+1)];
        Bv(i) = Kr*evalfr(Pi,x(i));
    else
        for j = 1:nNvs+1
            Av(i,j) = x(i)^(nNvs+1-j);
        end
        for j = 1:nNv+1
            Av(i,nNvs+1+j) = -evalfr(G1z,x(i))*x(i)^(nNv+1-j);
        end
    end
end
% Av = Av(3:end,:);
% Bv = Bv(3:end);
v = Av\Bv;
% v = pinv(Av)*Bv;
numVs = v(1:nNvs+1)';
numV  = v(nNvs+2:end)';
Vs = tf(numVs,denVs,Ts);
V  = tf(numV ,denV,Ts);

S   = phi-Vs*z^-d;
S   = minreal(S);
Ceq = 1/(S+1);
YR  = series(feedback(series(Ceq,Pr),V),Kr);
zpk(YR)
if (evalfr(V,1)-Kr < 1e-10)
    disp('V(1) == Kr, ok!')
else
    disp(['V(1) =/= Kr, error! ' num2str(evalfr(V,1)-Kr)])
end

eFlag = false;
if (evalfr(phi-Vs+1,1) < 1e-9)
    disp('Ceq(1) == inf, ok!')
else
    disp(['Ceq(1) =/= inf, error! ' num2str(evalfr(phi-Vs+1,1))])
%     eFlag = true;
end

[zz,pp,kk] = zpkdata(YR,'v');
for i = 1:Nrho
    if (abs(min(pp-rho(i))) < 1e-4)
        disp(['pole ' num2str(i) ' ok!'])
    else
        disp(['pole ' num2str(i) ' error! ' num2str(abs(min(pp-rho(i))))])
%         eFlag = true;
    end
end

disp(max(abs(pp)))
if(max(abs(pp))>=1)
    disp('There is a least one unstable pole!')
    eFlag = true;
end
% eFlag = false;
if ~eFlag
    % Call simulation
    sim('simuC1')
end