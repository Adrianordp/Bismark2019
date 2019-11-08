% mudar as incertezas no arquivo ex1_liu_geso

% clear all
close all
clc

Tsim = 22.5;
% ref = 0.2;
Truido = 20;
% Pruido = 0.000004;
% Aruido = 0.05;
Tpert = 10;
% Apert = -0.06;



s = tf('s');
Ts = 0.05;

G = 2/((3*s-1)*(s-1));
Lm = 0.3;
P = G*exp(-Lm*s)


Gp = G;
L = Lm;
% Gp = 2.1/((2.85*s-1)*(0.95*s-1));
% L = 0.315;
Pr = Gp*exp(-L*s)

Gpd = c2d(Gp,Ts);
d = round(L/Ts);

% dP = (Gp/G)*exp((Lm-L)*s)-1;


% d = round(L/Ts);
z = tf('z',Ts);
iz = inv(z);
Gd = c2d(G,Ts);
[nG,dG] =  tfdata(Gd,'v')
dn = round(Lm/Ts);
rd = roots(dG)


% lamb = 0.96

[A,C,B,D] = tf2ss(nG,dG)
sys = ss(A,C,B,D)
% csys = canon(sys,'companion')
A = A'
B = B'
C = C'

% Gd1 = ss(A,B,C,0,Ts);
% Gd2 = ss(A,B,eye(2),[0;0],Ts);

% A = [-dG(2:end)',[1;0]]
% B = nG(2:end)'
% C = [1 0];
% p1 = 0.96;
% p2 = 0.94;Gd1

% p1 = 0.9;
% p2 = 0.9;
% p1 = 0.95;
% p2 = 0.95;

% p1 = exp((-1/0.51)*Ts);
p1 = exp((-1/0.85)*Ts);
p2 = p1;

% mf = conv([1,-p1],[1,-p2])

% Mx = ctrb(A,B)
% PcA = A^2+mf(2)*A+mf(3)*eye(2);
% K = [0 1]*inv(Mx)*PcA
K = acker(A,B,[p1 p2])

%K = place(A,B,[p1 p2])
kr = inv([1 0]*inv(eye(2)-A+B*K)*B)
% kr = inv([1 0]*inv(B*K-A)*B)
% kr = inv([1 0]*inv(eye(2)-A+C*K)*C);


% g_ = K*inv(z*eye(2) - A)*B
% g_ = ss(A,B,K,0)
[ng_,dg_] = ss2tf(A,B,K,0)
g_ = tf(ng_,dg_,Ts)

ng_nG = tf(ng_,nG,Ts)

g1 = rd(1);
g2 = rd(2);

% Filtro1
% a1 = 0.94
% a2 = 0.94
% a3 = 0.7
a1 = 0.53
a2 = 0.53
a3 = 0.53
% d = 10;
af = [1 1 1;1 1/g1 1/g1^2;1 1/g2 1/g2^2;]
aux = ng_nG*(1-a1*iz)*(1-a2*iz)*(1-a3*iz)*z^dn;
% bf = [evalfr(aux,1);evalfr(aux,g1);evalfr(aux,g2)];
bf = [kr*(1-a1)*(1-a2)*(1-a3);evalfr(aux,g1);evalfr(aux,g2)];
x = inv(af)*bf
V = minreal((x(1)+x(2)*iz+x(3)*iz^2)/((1-a1*inv(z))*(1-a2*inv(z))*(1-a3*inv(z))));
[nV,dV] = tfdata(V,'v');

% 2 formas de se calcular S
% % Forma 1
% 
% GV = Gd*V;
% [nGV,dGV] = tfdata(GV,'v');
% [r,p,k] = residue(nGV,dGV)
% % Nv_ = deconv(nV,dG)
% Nv_ = r(3)*conv([1 -a2],[1 -a3])+r(4)*conv([1 -a1],[1 -a3])+r(5)*conv([1 -a1],[1 -a2]);
% 
% S_ = 0;
% for i = 1:d
%     S_ = S_ + K*A^(i-1)*B*iz^i;
% end
% S = minreal(S_ + tf(Nv_,dV,Ts))

% Forma 3
Snum = [conv(ng_,dV) zeros(1,dn)] - [zeros(1,dn) conv(nV,nG)]
Sden = [conv(dG,dV) zeros(1,dn)];
S = minreal(tf(Snum,Sden,Ts))

% % Forma 1 Bismark
% Ab = [];
% Bb = [];
% 
% xa = [g1,g2, a1, 0.4,0.7];  % dois ultimos valores randomicos
% for q = 1:5
%     x = xa(q);
%     Ab = [Ab;[[x 1]*polyval(dV,x),[x^2, x, 1]*polyval(dG,x)]];
%     Bb = [Bb;[polyval(nG,x)*polyval(nV,x)]];
% end
% 
% aux = Ab\Bb; 
% Nv_ = aux(3:end)'; 
% Ng_ = aux(1:2)';
% 
% aux1 = tf(Nv_,dV,Ts);
% aux2 = tf(Ng_,dG,Ts);  
% 
% S_ = 0;
% for ib = 1:(dn)
%     S_ =minreal(S_ + K*A^(ib-1)*B*iz^ib);
% end
% 
% S = minreal(parallel(S_,-aux1*iz^dn))


% ar1 = 0.95;
% ar2 = 0.9;

betaf1 = 0.8985;
betaf2 = betaf1;
alphaf1 = 0.94;
alphaf2 = alphaf1;

Frp = tf([(1-betaf1)*(1-betaf2),0,0],[1,-(betaf1+betaf2),betaf1*betaf2],Ts);
Frz = tf([1,-(alphaf1+alphaf2),alphaf1*alphaf2],[(1-alphaf1)*(1-alphaf2),0,0],Ts);
Fr = minreal(Frp*Frz);

sim('nonmininumS_2018a')

% robustez

X2 = log10(pi/Ts);
X1 = -2;
N = 100;


w = logspace(X1, X2, N);
jw = w*sqrt(-1);

zw= exp(jw*Ts);

for i = 1:length(zw)
    
    Pw = freqresp(P,w(i));
    Prw = freqresp(Pr,w(i));
    
    Gdw = freqresp(Gd,w(i));
    Vw = freqresp(V,w(i));
    %dPw(i) = abs(freqresp(dP,w(i))-1);
    dPw(i) = abs(Prw/Pw-1);
Ir_bis(i) = abs((1+K*inv(zw(i)*eye(length(A))-A)*B)/(Gdw*Vw));
end


fu

e = 1 - y;
efu = 1 - y_fu;

IAEsp = trapz(abs(e(1:Tpert/Ts)))*Ts
IAEspfu = trapz(abs(efu(1:Tpert/Ts)))*Ts

IAEld = trapz(abs(e(Tpert/Ts:Truido/Ts)))*Ts
IAEldfu = trapz(abs(efu(Tpert/Ts:Truido/Ts)))*Ts

TVld = sum(abs(diff(u(Tpert/Ts:Truido/Ts))))
TVldfu = sum(abs(diff(u_fu(Tpert/Ts:Truido/Ts))))

Var = var(u(Truido/Ts:end))
Varfu = var(u_fu(Truido/Ts:end))


f1 = figure;
subplot(2,1,1)
%plot(t,y,t,yliu,'r--','linewidth',2)
plot(t,y_fu,'linewidth',2)
hold on
plot(t,y,'r','linewidth',1)

axis([0 Tsim -0.2 1.2])

set(gca,'FontName','Times New Roman','FontSize',12)
leg1 = legend({'Fu and Tan [16]','Proposed'},'FontName','Times New Roman','FontSize',12,'location','best');
set(leg1(1),'Interpreter','latex');
grid
xlabel('Time (s)','interpreter','Latex')
ylabel('Output','interpreter','Latex')


subplot(2,1,2)

%plot(t,u,t,uliu,'r--','linewidth',2)
plot(t,u_fu,'linewidth',2)
hold on
plot(t,u,'r','linewidth',1)
axis([0 Tsim -5 10])

set(gca,'FontName','Times New Roman','FontSize',12)
%legend({'Proposed','Geng et al. [XX]'},'FontName','Times New Roman','FontSize',14,'location','best')
grid
xlabel('Time (s)','interpreter','Latex')
ylabel('Control Input','interpreter','Latex')

%Position: [403 246 560 420]
set(f1,'Position',[100 66 560 600])
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f1,'ex2_nom','-dpdf','-r0')
% print(f1,'ex2_pert','-dpdf','-r0')


f2 = figure(2);
% loglog(w,Ir_fu,'b','linewidth',2)
loglog(w,Ir_fu,'linewidth',2)
hold on
loglog(w,Ir_bis,'r','linewidth',1)

loglog(w,dPw,'k--','linewidth',2)
%loglog(w,Ir_bis,w,Ir_geng,w,dPw)


set(gca,'FontName','Times New Roman','FontSize',12);
%set(gca(1),'DisplayName','Proposed');

leg2 = legend({'Fu and Tan [16]','Proposed','$\overline{\delta P}$'},'FontName','Times New Roman','FontSize',12,'location','best');
set(leg2(1),'Interpreter','latex');




%legend({'Proposed','Geng et al. [XX]','$\overline{\delta P}'},'Interpreter','Latex','FontName','Times New Roman','FontSize',14,'location','best')

%legend({'$\overline{\delta P}$'},'Interpreter','latex');



grid
xlabel('Frequency $\omega$ (rad/s)','interpreter','Latex')
ylabel('Robustness Index $I_r$','interpreter','Latex')

axis([0.01 pi/Ts 0.1 100])

%Position: [403 246 560 420]
%set(f2,'Position',[700 246 560 420])

%set(f2,'Position',[403 246 560 420])

set(f2,'Units','Inches');
pos = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f2,'ex2_robustness','-dpdf','-r0')

