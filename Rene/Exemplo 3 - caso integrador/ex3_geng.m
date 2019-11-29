% mudar as incertezas no arquivo ex1_liu_geso

clear all
close all
clc

Tsim = 40;
ref = 0.2;
Truido = 35;
Pruido = 0.000004;
Aruido = 0.05;
Tpert = 20;
Apert = -0.06;

s = tf('s');
G = (1-0.2*s)/(s*(s-1));
L = 0.2;
Ts = 0.01;

wn = 20;
xi = 0.25;
% Pr = G*exp(-0.2*1.3*s)% Planta com incerteza
dP = exp(-0.06*s);
Pr = G*exp(-0.2*s);                                  % Sem incerteza

z = tf('z',Ts);
iz = inv(z);
Gd = c2d(G,Ts);
[nG,dG] =  tfdata(Gd,'v')

rd = roots(dG);
lamb = 0.93;        % Gang filter parameter

% State-space model
A = [-dG(2:end)',[1;0]]
B = nG(2:end)'
C = [1 0];

Gd1 = ss(A,B,C,0,Ts)
Gd2 = ss(A,B,eye(2),[0;0],Ts)
% 
% p1 = 0.95;
% p2 = 0.95;

p1 = 0.95;
p2 = 0.95;


a1 = 0.99;
a2 = a1;
a3 = a1;

K = acker(A,B,[p1 p2])

kr = inv([1 0]*inv(eye(2)-A+B*K)*B)

[ng_,dg_] = ss2tf(A,B,K,0)
g_ = tf(ng_,dg_,Ts)

ng_nG = tf(ng_,nG,Ts)

g1 = rd(1);
g2 = rd(2);

% Filtro1

% a1 = 0.988;
% a2 = 0.988;
% a3 = 0.99;


d = L/Ts;

syms w v0 v1 v2 numV

numG_ = poly2sym(ng_,w);
denG_ = poly2sym(dG,w);
numG = poly2sym(nG,w);
denG = poly2sym(dG,w);
numV = v0 + v1*w^-1 + v2*w^-2
denV = (1-a1*w^-1)*(1-a2*w^-1)*(1-a3*w^-1)

S_1 = denG*denV + numG_*denV - numG*numV*w^-d;
S_diff = diff(S_1,w)

num = subs(S_diff,w,1)

v0af = coeffs(num,v0);
v1af = coeffs(num,v1);
v2af = coeffs(num,v2);
coefbf = double(coeffs(num));

af = [1 1 1;1 1/g1 1/g1^2;double(v0af(end)) double(v1af(end)) double(v2af(end))]
aux = ng_nG*(1-a1*iz)*(1-a2*iz)*(1-a3*iz)*z^d;
bf = [evalfr(aux,1);
    evalfr(aux,g1);
    -coefbf(1)];
x = inv(af)*bf
V = minreal((x(1)+x(2)*iz+x(3)*iz^2)/((1-a1*inv(z))*(1-a2*inv(z))*(1-a3*inv(z))))
%bf
evalfr(V,1)

[nV,dV] = tfdata(V,'v');

% S_tfa = minreal(g_-Gd*V*z^-d);

nS = [conv(ng_,dV) zeros(1,d)] - [zeros(1,d) conv(nG,nV)];
dS = [conv(dG,dV) zeros(1,d)];
S_tf = minreal(tf(nS,dS,Ts));
% 3 formas de se calcular S
% Forma 1
% GV = Gd*V;
% [nGV,dGV] = tfdata(GV,'v');
% [r,p,k] = residue(nGV,dGV)
% % Nv_ = deconv(nV,dG)
% Nv_ = r(3)*conv([1 -a2],[1 -a3])+r(4)*conv([1 -a1],[1 -a3])+r(5)*conv([1 -a1],[1 -a2])
% S_ = 0;
% for i = 1:d
%     S_ = S_ + K*A^(i-1)*B*iz^i;
% end
%S = minreal(S_ + tf(Nv_,dV,Ts))
% 

% Forma 1 Bismark
Ab = [];
Bb = [];


xa = [g1,g2, a1, 0.4,0.7];  % dois ultimos valores randomicos
for q = 1:5
    x = xa(q);
    Ab = [Ab;[[x 1]*polyval(dV,x),[x^2, x, 1]*polyval(dG,x)]];
    Bb = [Bb;[polyval(nG,x)*polyval(nV,x)]];
end

aux = Ab\Bb; 
Nv_ = aux(3:end)'
Ng_ = aux(1:2)';

aux1 = tf(Nv_,dV,Ts);
aux2 = tf(Ng_,dG,Ts);  

S_ = 0;
for ib = 1:(d)
    S_ =minreal(S_ + K*A^(ib-1)*B*iz^ib);
end

S = minreal(parallel(S_,-aux1*iz^d))

% % Forma 2
% S = minreal(g_ - V*Gd*z^-d)
% 

% Forma 3
% Snum = [conv(ng_,dV) zeros(1,d)] - [zeros(1,d) conv(nV,nG)]
% Sden = [conv(dG,dV) zeros(1,d)];
% S = minreal(tf(Snum,Sden,Ts))

% Forma 4
% Internal instability solution
% S_til = (kc-Vd*z^-dn);
% [NS_til,DS_til] = tfdata(S_til,'v');
% T_til = deconv(NS_til,[1 -pnd]);
% S = tf(conv(Gndnum,T_til),DS_til,Ts);

% Internal instability solution
% T_til = deconv(Snum,dG);
% DS_til = [dV zeros(1,3)];
% S = tf(T_til(1:end-1),DS_til,Ts)

% ar1 = 0.95;
% ar2 = 0.9;

ar1 = 0.994;
ar2 = ar1;

Fr =  tf([(1-ar1)*(1-ar2),0,0],[1,-(ar1+ar2),ar1*ar2],Ts);

sim('sim1_nonmininumS')


%%
% robustez
X2 = log10(pi/Ts);
X1 = -2;
N = 100;

w = logspace(X1, X2, N);
jw = w*sqrt(-1);

zw= exp(jw*Ts);

for i = 1:length(zw)
    

    Gdw = freqresp(Gd,w(i));
    Vw = freqresp(V,w(i));
    dPw(i) = abs(freqresp(dP,w(i))-1);

Ir_bis(i) = abs((1+K*inv(zw(i)*eye(length(A))-A)*B)/(Gdw*Vw));
end

ex3_liu_geso

f1 = figure;
subplot(2,1,1)
%plot(t,y,t,yliu,'r--','linewidth',2)
plot(t,yliu,'linewidth',2)
hold on
plot(t,y,'r','linewidth',1)

axis([0 Tsim -0.05 0.25])

set(gca,'FontName','Times New Roman','FontSize',12)
leg1 = legend({'Geng et al. [17]','Proposed'},'FontName','Times New Roman','FontSize',12,'location','best');
set(leg1(1),'Interpreter','latex');
grid
xlabel('Time (s)','Interpreter','latex')
ylabel('Output','Interpreter','latex')



subplot(2,1,2)

%plot(t,u,t,uliu,'r--','linewidth',2)
plot(t,uliu,'linewidth',2)
hold on
plot(t,u,'r','linewidth',1)
axis([0 Tsim -0.15 0.25])

set(gca,'FontName','Times New Roman','FontSize',12)
%legend({'Proposed','Geng et al. [XX]'},'FontName','Times New Roman','FontSize',14,'location','best')
grid
xlabel('Time (s)','Interpreter','latex')
ylabel('Control Input','Interpreter','latex')

%Position: [403 246 560 420]
set(f1,'Position',[100 66 560 600])
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f1,'ex3_nom','-dpdf','-r0')
% print(f1,'ex3_perturbed','-dpdf','-r0')

f2 = figure(2);
loglog(w,Ir_geng,'linewidth',2)
hold on
loglog(w,Ir_bis,'r','linewidth',1)

loglog(w,dPw,'k--','linewidth',2)
%loglog(w,Ir_bis,w,Ir_geng,w,dPw)


set(gca,'FontName','Times New Roman','FontSize',12)
%set(gca(1),'DisplayName','Proposed');

leg2 = legend({'Geng et al. [17]','Proposed','$\overline{\delta P}$'},'FontName','Times New Roman','FontSize',12,'location','best');
set(leg2(1),'Interpreter','latex');

%legend({'Proposed','Geng et al. [XX]','$\overline{\delta P}'},'Interpreter','Latex','FontName','Times New Roman','FontSize',14,'location','best')

%legend({'$\overline{\delta P}$'},'Interpreter','latex');

grid
xlabel('Frequency $\omega$ (rad/s)','Interpreter','latex')
ylabel('Robustness Index $I_r$','interpreter','Latex')

axis([0.01 pi/Ts 0.1 100])
set(f2,'Units','Inches');
pos = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f2,'ex3_robustness','-dpdf','-r0')

% Indices calculation
e = 0.2 - y;
eliu = 0.2 - yliu;

IAE = trapz(abs(e(Tpert/Ts:Truido/Ts)))*Ts
IAEgeng = trapz(abs(eliu(Tpert/Ts:Truido/Ts)))*Ts

TV = sum(abs(diff(u(Tpert/Ts:Truido/Ts))))
TVgeng = sum(abs(diff(uliu(Tpert/Ts:Truido/Ts))))

Var = var(u(Truido/Ts:end))
Vargeng = var(uliu(Truido/Ts:end))
