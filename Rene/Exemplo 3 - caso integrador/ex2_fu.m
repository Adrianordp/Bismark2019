% mudar as incertezas no arquivo ex1_liu_geso

clear all
close all
clc

s = tf('s');
G = 2/((3*s-1)*(s-1));

% G = tf([-1,1],[6 5 1]);
% G = zpk(G);
Ts = 0.05;

wn = 20;
xi = 0.25;
% Pr = G*exp(-1.25*s)*(wn^2/(s^2+2*xi*s*wn+wn^2)) % Planta com incerteza
Pr = G*exp(-0.3*s);                                  % Sem incerteza

z = tf('z',Ts);
iz = inv(z);
Gd = c2d(G,Ts);
[nG,dG] =  tfdata(Gd,'v')

rd = roots(dG)


% lamb = 0.96

% [A,B,C,D] = tf2ss(nG,dG)
% sys = ss(A,B,C,D);
% csys = canon(sys,'companion')
% A = A'
% B = B'
% C = C'
A = [-dG(2:end)',[1;0]]
B = nG(2:end)'
C = [1 0];
% p1 = 0.96;
% p2 = 0.94;

% p1 = 0.9;
% p2 = 0.9;
p1 = 0.95;
p2 = 0.95;

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
a1 = 0.97
a2 = 0.97
a3 = 0.97
% d = 10;
d = 6;
af = [1 1 1;1 1/g1 1/g1^2;1 1/g2 1/g2^2;]
aux = ng_nG*(1-a1*iz)*(1-a2*iz)*(1-a3*iz)*z^d;
bf = [kr*(1-a1)*(1-a2)*(1-a3);evalfr(aux,g1);evalfr(aux,g2)];
x = inv(af)*bf
V = minreal((x(1)+x(2)*iz+x(3)*iz^2)/((1-a1*inv(z))*(1-a2*inv(z))*(1-a3*inv(z))));
%bf
evalfr(V,1)

% 3 formas de se calcular S
% Forma 1
[nV,dV] = tfdata(V,'v');

GV = Gd*V;
[nGV,dGV] = tfdata(GV,'v');
[r,p,k] = residue(nGV,dGV)
% Nv_ = deconv(nV,dG)
Nv_ = r(3)*conv([1 -a2],[1 -a3])+r(4)*conv([1 -a1],[1 -a3])+r(5)*conv([1 -a1],[1 -a2]);

% S_ = 0;
% for i = 1:d
%     S_ = S_ + K*A^(i-1)*B*iz^i;
% end
%S = minreal(S_ + tf(Nv_,dV,Ts))
% 
% % Forma 1 Bismark
% 
% Ab = [];
% Bb = [];
% 
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
% for ib = 1:(d)
%     S_ =minreal(S_ + K*A^(ib-1)*B*iz^ib);
% end
% 
% S = minreal(parallel(S_,-aux1*iz^d))


% % Forma 2
% S = minreal(g_ - V*Gd*z^-d)
% 
% % Forma 3
Snum = [conv(ng_,dV) zeros(1,d)] - [zeros(1,d) conv(nV,nG)]
Sden = [conv(dG,dV) zeros(1,d)];
S = minreal(tf(Snum,Sden,Ts))

% ar1 = 0.95;
% ar2 = 0.9;

ar1 = 0.97
ar2 = ar1;

Fr =  tf([(1-ar1)*(1-ar2),0,0],[1,-(ar1+ar2),ar1*ar2],Ts);

sim('sim1_nonmininumS')

% ex1_liu_geso

% e = 1 - y;
% % eliu = 1 - yliu;
% 
% IAE = trapz(abs(e(300:550)))*Ts
% % IAEgeng = trapz(abs(eliu(300:550)))*Ts
% 
% TV = sum(abs(diff(u(300:550))))
% % TVgeng = sum(abs(diff(uliu(300:550))))
% 
% Var = var(u(550:end))
% % Vargeng = var(uliu(550:end))

subplot(2,1,1)
%plot(t,y,t,yliu,'r--','linewidth',2)
plot(t,y,'linewidth',2)
% hold on
% plot(t,yliu,'r','linewidth',1)

%axis([0 60 -0.2 1.2])

% set(gca,'FontName','Times New Roman','FontSize',14)
% legend({'Proposed','Geng et al. [XX]'},'FontName','Times New Roman','FontSize',14,'location','best')
% grid
xlabel('Time (s)')
ylabel('Output')


subplot(2,1,2)

%plot(t,u,t,uliu,'r--','linewidth',2)
plot(t,u,'linewidth',2)
% hold on
% plot(t,uliu,'r','linewidth',1)
axis([0 60 -0 2.5])

% set(gca,'FontName','Times New Roman','FontSize',14)
% %legend({'Proposed','Geng et al. [XX]'},'FontName','Times New Roman','FontSize',14,'location','best')
% grid
xlabel('Time (s)')
ylabel('Control Input')

