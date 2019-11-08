% mudar as incertezas no arquivo ex1_liu_geso

clear all
close all
clc

s = tf('s');
%G = (-s+1)/((3*s+1)*(2*s+1));

G = tf(2,[3 -4 1]);
G = zpk(G);
Ts = 0.1;

wn = 20;
xi = 0.25;
Pr = G*exp(-1.25*s)*(wn^2/(s^2+2*xi*s*wn+wn^2)) % Planta com incerteza
Pr = G*exp(-s);                                  % Sem incerteza

dP = exp(-0.25*s)*(wn^2/(s^2+2*xi*s*wn+wn^2))%       incerteza multiplicativa

z = tf('z',Ts);
iz = inv(z);
Gd = c2d(G,Ts);
[nG,dG] =  tfdata(Gd,'v')

rd = roots(dG);


lamb = 0.96


A = [-dG(2:end)',[1;0]];
B = nG(2:end)';
C = [1 0];
% p1 = 0.96;
% p2 = 0.94;

p1 = 0.9;
p2 = 0.9;

mf = conv([1,-p1],[1,-p2])

Mx = ctrb(A,B)
PcA = A^2+mf(2)*A+mf(3)*eye(2);
K = [0 1]*inv(Mx)*PcA

%K = place(A,B,[p1 p2])
kr = inv(C*inv(eye(2)-A+B*K)*B);


g1 = rd(1);
g2 = rd(2);

% Filtro1
% a1 = 0.94
% a2 = 0.94
% a3 = 0.7

a1 = 0.93
a2 = 0.93
a3 = 0.8
d = 10;
af = [1 1 1;1 1/g1 1/g1^2;1 1/g2 1/g2^2;]
%aux = (1-a1*iz)*((1-a2*iz)*(1-a3*iz))*z^d;
%bf = [kr*(1-a1)*(1-a2)*(1-a3);evalfr(aux,p1);evalfr(aux,p1)];  % versão
%anterior
Ga = ss(A,B,K,0,Ts);
[nGa, dGa] = tfdata(Ga,'v')

aux2 = tf(nGa,nG,Ts)
aux = aux2*(1-a1*iz)*((1-a2*iz)*(1-a3*iz))*z^d;
bf = [kr*(1-a1)*(1-a2)*(1-a3);evalfr(aux,g1);evalfr(aux,g1)];
x = inv(af)*bf
V = minreal((x(1)+x(2)*iz+x(3)*iz^2)/((1-a1*inv(z))*(1-a2*inv(z))*(1-a3*inv(z))));
%bf
evalfr(V,1)

% ar1 = 0.95;
% ar2 = 0.9;

ar1 = 0.91
ar2 = ar1;


Fr =  tf([(1-ar1)*(1-ar2),0,0],[1,-(ar1+ar2),ar1*ar2],Ts);
S = Ga-Gd*V*z^(-d);
sim('sim2_nonmininum.slx')
plot(t,y,'linewidth',2)

axis([0 60 -0.2 1.2])

set(gca,'FontName','Times New Roman','FontSize',14)
% legend({'Proposed','Geng et al. [XX]'},'FontName','Times New Roman','FontSize',14,'location','best')
grid
xlabel('Time (s)')
ylabel('Output')

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

ex1_liu_geso


subplot(2,1,1)
%plot(t,y,t,yliu,'r--','linewidth',2)
plot(t,y,'linewidth',2)
hold on
plot(t,yliu,'r','linewidth',1)

axis([0 60 -0.2 1.2])

set(gca,'FontName','Times New Roman','FontSize',14)
legend({'Proposed','Geng et al. [XX]'},'FontName','Times New Roman','FontSize',14,'location','best')
grid
xlabel('Time (s)')
ylabel('Output')



subplot(2,1,2)

%plot(t,u,t,uliu,'r--','linewidth',2)
plot(t,u,'linewidth',2)
hold on
plot(t,uliu,'r','linewidth',1)
axis([0 60 -0 2.5])

set(gca,'FontName','Times New Roman','FontSize',14)
%legend({'Proposed','Geng et al. [XX]'},'FontName','Times New Roman','FontSize',14,'location','best')
grid
xlabel('Time (s)')
ylabel('Control Input')


figure
loglog(w,Ir_bis,'linewidth',2)
hold on
loglog(w,Ir_geng,'r','linewidth',1)

loglog(w,dPw,'k--','linewidth',2)
%loglog(w,Ir_bis,w,Ir_geng,w,dPw)


set(gca,'FontName','Times New Roman','FontSize',14)
%set(gca(1),'DisplayName','Proposed');

leg1 = legend({'Proposed','Geng et al. [XX]','$\overline{\delta P}$'},'FontName','Times New Roman','FontSize',14,'location','best')
set(leg1(1),'Interpreter','latex');


%legend({'Proposed','Geng et al. [XX]','$\overline{\delta P}'},'Interpreter','Latex','FontName','Times New Roman','FontSize',14,'location','best')

%legend({'$\overline{\delta P}$'},'Interpreter','latex');



grid
xlabel('Frequency \omega (rad/s)')
ylabel('Robustness Index $I_r$','interpreter','Latex')


axis([0.01 pi/Ts 0.1 100])