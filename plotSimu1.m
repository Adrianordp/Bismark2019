% Author: Adriano Rodrigues
% CloseFcn callback for Example_2_1.slx and plot all relevant signals

% Close figures
close all

% Plot y(t)
hy = figure;               % Initialize figure frame
plot(t,r,'LineWidth',2)    % Plot curves
grid on                    % Enable grid
hold on
plot(t,y,'LineWidth',2)
axis([min(t) max(t) -.2 1.1]); % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('y(t)')
title('y(t)')
legend('r(t)','y(t)')
hy.Position = [-1066 506 895 385];

% Plot u(t)
hu = figure;              % Initialize figure frame
plot(t,u,'LineWidth',2)   % Plot curves
grid on                   % Enable grid
axis([min(t) max(t) 2.6]) % Limit axes
% Name labels
xlabel('Time [sample]')
ylabel('u(t)')
title('u(t)')
legend('u(t)')
hu.Position = [-1066 506 895 385];