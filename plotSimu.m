% Author: Adriano Rodrigues
% CloseFcn callback for Example_2_1.slx and plot all relevant signals

% Close figures
close all

% Plot y(t)
figure                     % Initialize figure frame
plot(t,y,'LineWidth',2)    % Plot curves
grid on                    % Enable grid
hold on
plot(t,r,'LineWidth',2)
% axis([min(t) max(t) -3 1]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('y(t)')
title('y(t)')
legend('y(t)','r(t)')

% Plot u(t)
figure                      % Initialize figure frame
plot(t,u,'LineWidth',2)     % Plot curves
grid on                     % Enable grid
% axis([min(t) max(t) -1 .2]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('u(t)')
title('u(t)')
legend('u(t)')