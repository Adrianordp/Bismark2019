% descomentar at� lamb
% % clear all
% % close all
% % clc
% % 
% % s = tf('s');
% % G = tf([-1,1],[6 5 1]);
% % G = zpk(G);
% % Ts = 0.1;
% % 
% % 
% % wn = 20;
% % xi = 0.25;
% % 
% % Pr = G*exp(-1.25*s)*(wn^2/(s^2+2*xi*s*wn+wn^2));
% % %Pr = G*exp(-s);
% % 
% % 
% % 
% % 
% % z = tf('z',Ts);
% % iz = inv(z);
% % Gd = c2d(G,Ts);
% % [nG,dG] =  tfdata(Gd,'v');
% % 
% % rd = roots(dG);
% % 
% % 
% % lamb = 0.96;

% eq 13
nk = 0;              % for noise rejection
Hl = ((1-lamb)*z)^nk*(z-1)/(z-lamb)^(nk+1);

% eq 12
m = 1;

aux1 = ((z-lamb)^m/(dG(1)*z^2+dG(2)*z+dG(3)));
aux2 = minreal(inv(Hl));
Gt = minreal(((z-lamb)^m/(dG(1)*z^2+dG(2)*z+dG(3)))*inv(Hl));
[nt,dt] = tfdata(Gt,'v');

 Am = [0 1  ;
       -dG(end:-1:2)];
  Bm = [0;nG(end)];
  Cm = [1,nG(2)/nG(3)];
  
  Ae = [Am,[0;1];[0 0 1]];
  Be = [Bm;0];
  Ee = [0;0;1];
  
 Ce = [Cm,0];
  
 b0 = 0.9815;
 bc = 0.99;              % No artigo este valor foi arredondado

 aux = conv([1 -2*b0 b0^2],[1,-b0]);
 fia = Ae^3 +  aux(2)*Ae^2 + aux(3)*Ae +  aux(4)*eye(3);
 
 L0  = fia*inv([Ce;Ce*Ae;Ce*Ae^2])*[0;0;1];
 
 % ?0.9065 0.9530
 k1 = (-bc)^2-dG(3);
 k2 = 2*(-bc)^1-dG(2);
 
 K0b = [k1,k2]/nG(end);
 K0 = [K0b,1/nG(end)];
 
% lf nf 0.94 2 
% 0.2146(z ? 0.9512)2z
% (z ? 0.9047)(z ? 0.94)2
lf = 0.992;
nf = 3;

aux = zpk(ss(Am-Bm*K0b,Bm,Cm,0,Ts));

% Cm*inv(z*eye(2)-Am+Bm*K0b)*Bm

%Tm = ()

Kf = 0.00024875*(z-bc)^2/((z-0.9512)*(z-lf)^nf);

Ap = [ -dt(2:end);1 0 0; 0 1 0];
Bp = [1;0;0];
Cp = [nt(2:end)];
 
 % eq 14
 d0 = d;
 
 Gta = minreal(ss(Ap,Ap^d0*Bp,Cp,0,Ts));
 
 [ngta,dgta] = tfdata(Gta,'v');
 
 F2 = tf(ngta,conv([1,-lamb],[1,-lamb]),Ts);
 
 % Fir filter
 clear aux
 for k = 1:d0
    aux(k) = Cp*Ap^(k-1)*Bp;
 end

 F1 = tf(aux,[1,zeros(1,d0)],Ts);
 Tau = tf(nG,[1,-lamb],Ts);
 Tau = Tau*Hl;
 
 F1 = F1*Tau;

A = [-dG(2:end)',[1;0]];
B = nG(2:end)';
% p1 = 0.96;
% p2 = 0.94;

p1 = 0.9;
p2 = 0.85;

K = place(A,B,[p1 p2]);
kr = inv([1 0]*inv(eye(2)-A+B*K)*B);

sim('sim3_nonmininum_geso')


% % % robustez
% % 
% X2 = log10(pi/Ts);
% X1 = -2;
% N = 100;

w = logspace(X1, X2, N);
jw = w*sqrt(-1);

zw= exp(jw*Ts);

for i = 1:length(zw)
    
    F1w = freqresp(F1,w(i));
    F2w = freqresp(F2,w(i));
    Pdw = freqresp(Gd*z^-d0,w(i));

aux1 =1+ K0*inv(zw(i)*eye(length(Ae))+L0*Ce-Ae)*(Be+L0*F1w);
aux2 = K0*inv(zw(i)*eye(length(Ae))+L0*Ce-Ae)*L0*F2w;
aux3 = aux2/aux1;

Ir_geng(i) = abs((1+Pdw*aux3)/(Pdw*aux3));
end
%  figure
%  loglog(w,Ir_geng)
