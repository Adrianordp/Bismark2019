reseta

Ts = 0.1;
s = tf('s');
z = tf('z',Ts);
P1s = (-s+1)/(2*s+1)/(3*s+1)*exp(-s);
step(P1s)

P1z = c2d(P1s,Ts);
d = P1z.InputDelay;

K = [39.5284 40.4371];
V = (5.155*z^3 - 9.848*z^2 + 4.704*z)/(z-0.9)^2/(z-0.74);

[A,B,C,D]  = tf2ss(P1z.num{1},P1z.den{1});
A = A';
buff = B';
B = C';
C = buff;

beta = [0.92 0.92];

K = acker(A,B,beta)
Kr = inv(C/(eye(size(A))+B*K-A)*B);
F = Kr*(1-0.89)^2*z^2/(z-0.89)^2;

sim('Bismark2019Simu')

