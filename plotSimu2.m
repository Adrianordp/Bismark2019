% Author: Adriano Rodrigues
% CloseFcn callback for Example_2_1.slx and plot all relevant signals

% Close figures
close all

% Plot y(t)
figure                     % Initialize figure frame
plot(t,r,'LineWidth',2)    % Plot curves
grid on                    % Enable grid
hold on
plot(t,y,'LineWidth',2)
axis([min(t) max(t) -.2 1.2]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('y(t)')
title('y(t)')
legend('r(t)','y(t)')
hy = gcf;
hy.Position = [-1056 581 875 395];

% Plot u(t)
figure                      % Initialize figure frame
plot(t,u,'LineWidth',2)     % Plot curves
grid on                     % Enable grid
axis([min(t) max(t) -5 10]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('u(t)')
title('u(t)')
legend('u(t)')
hu = gcf;
hu.Position = [-1056 581 875 395]