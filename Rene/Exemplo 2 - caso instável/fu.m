
% Ts = 0.01;
% s = tf('s');

% G = tf(2,conv([3 -1],[1 -1]));
% G = 2/((3*s - 1)*(s - 1));
% Lm = 0.3;
[num,den] = pade(Lm,1);
delay = tf(num,den);
Pn = G*delay;
[numG,denG] = tfdata(Pn,'v');

% L = 0.3;
% Pr = G;
% Pr.outputdelay = L;

% numG = numG/denG(1);
% denG = denG/denG(1);
[A,B,C,D] = tf2ss(numG,denG);
% [A,B,C,D] = tf2ss(numG,denG);
% [A,C,B,D] = tf2ss(numG/denG(1),denG/denG(1));
% A = A';
% B = B';
% C = C';

A = fliplr(flipud(A));
% A = rot90(A,2)
% % B = fliplr(flipud(B));
% B = flipud(B)*numG(end);
% C = fliplr(C)/numG(end);
b0 = C(end);
B = flipud(B)*b0;
C = fliplr(C)/b0;
% C = C';

N = 1/((0.51*s + 1)^2);
% N.outputdelay = Lm;
Cr = 0.5*((3*s - 1)*(s - 1))/((0.51*s + 1)^2);
K0 = [61341.2837 2914.1104 92.2411 1];
K0 = [acker(A,B,[-30 -30 -30]) 1];
% K0 = [place(A,B,[-30.00001 -30.00000001 -30]) 1];
% K0 = K0(1:end-1);
L0 = [0.4256 4.9673 2.9762 1.5846]';
% L0 = L0(1:end-1);
Bd = B*K0(end);
% Bd = B;

Ae = [A Bd; zeros(1,length(A)+1)];
Be = [B;0];
Ce = [C 0];
% Ae = [Bd A; zeros(1,length(A)+1)];
% Be = [B;0];
% Ae = A;
% Be = B;

% A = [0 1 0; 0 0 1; 0 0 0];
% B = [0 1 0]';
% 
% Ae = [A B; zeros(1,length(A)+1)];
% Be = [B;0];
% Ce = [C 0];
L0 = acker(Ae',Ce',[-2 -2 -2 -2])';
% L0 = place(Ae',Ce',[-2 -2 -2 -2])';

K = K0*inv(s*eye(4) - Ae + Be*K0 + L0*Ce)*L0;

sim fu_adrc2_2018a

% robustez

X2 = log10(pi/Ts);
X1 = -2;
N = 100;


w = logspace(X1, X2, N);
jw = w*sqrt(-1);

% zw= exp(jw*Ts);

% for i = 1:length(zw)
for i = 1:length(jw)
    

    Pw = freqresp(P,w(i));
    Kw = freqresp(K,w(i));
%     dPw(i) = abs(freqresp(dP,w(i))-1);

Ir_fu(i) = abs((1+Kw*Pw)/(Kw*Pw));
end




