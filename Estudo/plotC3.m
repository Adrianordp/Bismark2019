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
axis([min(t) max(t) -.05 0.25]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('y(t)')
title('y(t)')
legend('r(t)','y(t)')
hy = gcf;
hy.Position = [-1065 524 892 396];

% Plot u(t)
figure                      % Initialize figure frame
plot(t,u,'LineWidth',2)     % Plot curves
grid on                     % Enable grid
axis([min(t) max(t) 0 .3]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('u(t)')
title('u(t)')
legend('u(t)')
hu = gcf;
hu.Position = [-1065 524 892 396];