% Author: Adriano Rodrigues
% CloseFcn callback for Example_2_1.slx and plot all relevant signals

% Close figures
close all

% Plot u(t)
figure                      % Initialize figure frame
plot(t,u,'LineWidth',2)     % Plot curves
axis([min(t) max(t) 0 2.6]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('u(t)')
title('u(t)')
legend('u(t)')
hu = gcf;
hu.Position = [684 1 683 661];
grid on                     % Enable grid

% Plot y(t)
figure                     % Initialize figure frame
plot(t,r,'LineWidth',2)    % Plot curves
hold on
plot(t,y,'LineWidth',2)
% axis([min(t) max(t) -.05 1.2]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('y(t)')
title('y(t)')
legend('r(t)','y(t)')
hy = gcf;
hy.Position = [1 1 683 661];
step(P1s,Tsim)
axis([min(t) max(t) -.2 1.1]) % Limit axes
grid on                    % Enable grid