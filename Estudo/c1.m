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
r1     = .92; % .975
r2     = .92; % .97
rho    = [r1 r2];
alphaf = .00; % .95
betaf  = .89; % .89
betav  = [.9 .9 .74]; % .95
betavs = [.9 .9 .74]; % .965
Nrho   = length(rho);
Np     = length(p);
nz     = 2;
nDv    = length(betav);
nDvs   = length(betavs);
nNv    = 3;
nNvs   = 1;

% Control gain
K = acker(A,B,rho);

% Reference filter
Kr = 1/(C/(eye(size(A))+B*K-A)*B);
% Kr = 4.0028
F  = minreal(Kr*(1-betaf)^2*z^2/(z-betaf)^2*(1-alphaf*z^-1)^nz/(1-alphaf)^nz);
% F = 
%        0.04843 z^2
%   ---------------------
%   z^2 - 1.78 z + 0.7921

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

GaG = tf(nGa,nG,Ts);

phi = 0;
for i = 1:d
    phi = phi + K*A^(i-1)*B*z^-i;
end
phi = minreal(phi);
Nr  = nNvs+nNv-Nrho-Np;

x  = [1 p 0.8114];
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
    if i == 1
        for j = 1:nNv+1
            Av(i,j+nNvs+1) = x(i)^(nNv-j+1);
        end
        Bv(i) = Kr*evalfr(Pi,x(i));
    elseif i == Nx
        for j = 1:nNv+1
            Av(i,j+nNvs+1) = x(i)^(nNv-j+1);
        end
        Bv(i) = Kr*evalfr(Pi,x(i));
    else
        for j = 1:nNv+1
            Av(i,j+nNvs+1) = x(i)^(nNv-j+1);
        end
        Bv(i) = evalfr(Pi,x(i))*evalfr(GaG,x(i))*x(i)^d;
    end
end
Avb = Av;
Av = Av(1:3,3:5);
% Av = Av(:,3:6);
% Av(4,:) = [1 -1 1 -1];
Bv = Bv(1:end-1)
% Bv = Bv(1:end);
% Bv(4) = 0;
v = Av\Bv;
numV  = [v' 0];
% numV  = v';
% Vs = tf(numVs,denVs,Ts);
V  = tf(numV,denV,Ts);
%        5.155 z^2 - 9.848 z + 4.704
% V = ---------------------------------
%     z^3 - 2.54 z^2 + 2.142 z - 0.5994

x  = rand(1,5);
% x = rand(1,4);
% x = [beta1 beta2 beta3];
Nx = length(x);
Aa = zeros(Nx);
b  = zeros(Nx,1);
for i = 1:Nx
    dGZ = polyval(dG,x(i));
    dVZ = polyval(denV,x(i));
    nGZ = polyval(nG,x(i));
    nVZ = polyval(numV,x(i));
    Aa(i,:) = [[x(i)^2 x(i) 1]*dGZ [x(i) 1]*dVZ];
%     Aa(i,:) = [x(i)^2 x(i) 1]*dGZ;
    b(i) = nGZ*nVZ;
end

vs = Aa\b;

Vs = tf([vs(1) vs(2) vs(3)],denV,Ts);

S   = phi-Vs*z^-d;
S   = minreal(S);
Ceq = 1/(S+1);
YR  = series(feedback(series(Ceq,Pr),V),Kr);
zpk(YR)
YQ = Pr*(phi-Vs*z^-d+1)/(Pr*V+phi-Vs*z^-d+1);
if (evalfr(V,1)-Kr < 1e-10)
    disp('V(1) == Kr, ok!')
else
    disp(['V(1) =/= Kr, error! ' num2str(evalfr(V,1)-Kr)])
end

eFlag = false;
if (evalfr(S+1,1) < 1e-9)
    disp('Ceq(1) == inf, ok!')
else
    disp(['Ceq(1) =/= inf, error! ' num2str(evalfr(S+1,1))])
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
% if ~eFlag
%     % Call simulation
%     sim('simuC1')
% end

clear j
w  = logspace(-2,1.5,100);
Nw = length(w);
P  = zeros(Nw,1);
for i = 1:Nw
    ejw = exp(j*w(i)*Ts);
    P(i) = abs(1+K*((ejw*eye(size(A))-A)\B))/abs(freqresp(V,w(i))*freqresp(G1z,w(i)));
end
figure, loglog(w,P)
grid on

IR = gcf;
IR.Children.YLim = [.10 100];
IR.Children.XLim = [.01 10^1.5];
IR.Position = [-930 495 576 429];















