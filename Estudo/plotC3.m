% Author: Adriano Rodrigues
% CloseFcn callback for Example_2_1.slx and plot all relevant signals

% Close figures
close all

% Plot u(t)
figure                      % Initialize figure frame
hold on
plot(t,u,'r','LineWidth',2)     % Plot curves
plot(t,uu1,'--k','LineWidth',2)     % Plot curves
axis([min(t) max(t) 0 2.6]) % Limit axes
% Name labels
xlabel('Tempo [amostras]')
ylabel('u(t)')
title('u(t)')
legend('estudado','Torrico')
hu = gcf;
hu.Position = [-1122 388 1015 462];
grid on                     % Enable grid
hu.Children(2).GridAlpha = .3;
hu.Children(2).GridColor = [0 0 0];

% Plot y(t)
figure                     % Initialize figure frame
plot(t,r,'LineWidth',2)    % Plot curves
hold on
plot(t,y,'LineWidth',2)
plot(t,yy1,'--k','LineWidth',2)
% axis([min(t) max(t) -.05 1.2]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('y(t)')
title('y(t)')
legend('r(t)','estudado','Torrico')
hy = gcf;
hy.Position = [-1122 388 1015 462];
% step(P1s,Tsim)
axis([min(t) max(t) -.2 1.1]) % Limit axes
grid on                    % Enable grid
% min(y(200:end))
hy.Children(2).GridAlpha = .3;
hy.Children(2).GridColor = [0 0 0];