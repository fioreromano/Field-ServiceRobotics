%% 4a

clc
close all
clear all

% Parameters
g = 9.81;
l = 1.0;
alpha = pi/8;
gamma = 0.08;

w1 = sqrt(2*g/l*(1-cos(gamma-alpha)));
% Vettore di velocità angolari iniziali
thetadot0_vec = [-2.2, -1.4, -0.5, -0.1, 0.05, 0.5, 0.95, 1, 1.5, 10];
colors = lines(length(thetadot0_vec));

figure(1); hold on; grid on;
title('Phase Portraits for different $\dot{\theta}_0$', 'Interpreter', 'latex');
xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');


figure(2)
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('State', 'Interpreter', 'latex');
title('Rimless Wheel Dynamics over Time', 'Interpreter', 'latex');
legend show;
set(legend, 'Interpreter', 'latex');


figure(3); clf; hold on; grid on;
title('Time intervals between impacts for different $\dot{\theta}_0$', 'Interpreter', 'latex');
xlabel('Impact number', 'Interpreter', 'latex');
ylabel('Time interval $\Delta t$ (s)', 'Interpreter', 'latex');


for k = 1:length(thetadot0_vec)

    thetadot0 = thetadot0_vec(k);
    if thetadot0 >= 0
        theta0 = gamma - alpha;
    else
        theta0 = gamma + alpha;
    end

    y0 = [theta0; thetadot0];
    double_support = 0;

    T = [];
    Y = [];
    impactTimes = []; 
    t0 = 0; tf = 25; dt = 0.01;

    while t0 < tf
        options = odeset('Events', @(t, y) impact_event(t, y, alpha, gamma), 'MaxStep', dt);
        [t, y, te, ye, ie] = ode45(@(t, y) dynamics(t, y, g, l, double_support), [t0 tf], y0, options);

        T = [T; t];
        Y = [Y; y];

        if ~isempty(te)
            impactTimes = [impactTimes; te(end)]; 

            [y0, double_support] = impact_map(ye(end,:)', alpha, g, l); % ye(end,:) trasposto colonna
            t0 = te(end);   % usa l'ultimo impatto
        else
            break;
        end
    end

    % === PLOT ===
    figure(1)
    subplot(5, 2, k)
    plot(Y(:,1), Y(:,2), 'Color', colors(k,:), 'LineWidth', 1);
    hold on;
    plot(Y(1,1), Y(1,2), '*', 'Color', colors(k,:), 'LineWidth', 1.5);
    grid on;
    xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
    ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');
    title(['Phase Portrait for $\dot{\theta}_0 = $' num2str(thetadot0)], 'Interpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');   
    print(gcf, 'rimless_limit.pdf', '-dpdf', '-bestfit');
    
    figure(2)
    subplot(5, 2, k)
    plot(T, Y(:,1), '-', 'Color', colors(k,:), 'LineWidth', 0.8);
    hold on;
    plot(T, Y(:,2), '--', 'Color', colors(k,:), 'LineWidth', 0.8);
    grid on;
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('State', 'Interpreter', 'latex');
    title(['Dynamics for $\dot{\theta}_0 = $' num2str(thetadot0)], 'Interpreter', 'latex');
    legend({'$\theta$', '$\dot{\theta}$'}, 'Interpreter', 'latex', 'Location', 'best'); 
    set(gcf, 'PaperPositionMode', 'auto');   
    print(gcf, 'rimless_dyn.pdf', '-dpdf', '-bestfit');

    figure(3)
    if length(impactTimes) > 1
        delta_t = diff(impactTimes);
        plot(1:length(delta_t), delta_t, '-o', 'Color', colors(k,:), 'DisplayName', ['$\dot{\theta}_0=$' num2str(thetadot0)], 'LineWidth', 1);
        legend('Interpreter', 'latex', 'Location', 'best');
        set(gcf, 'PaperPositionMode', 'auto');   
        print(gcf, 'rimless_impact.pdf', '-dpdf', '-bestfit');
    else
        plot(1, NaN, 'x', 'Color', colors(k,:), 'DisplayName', ['$\dot{\theta}_0=$' num2str(thetadot0)]);
        legend('Interpreter', 'latex', 'Location', 'best');
        set(gcf, 'PaperPositionMode', 'auto');   
        print(gcf, 'rimless_impact.pdf', '-dpdf', '-bestfit');
    end

end

%% 
clc; close all; clear all;

% Parametri
g = 9.81;
l = 1.0;
alpha = pi/8;
gamma = 0.08;

thetadot0_vec = linspace(-10, 10, 500);
colors = zeros(length(thetadot0_vec), 3); 

for k = 1:length(thetadot0_vec)
    thetadot0 = thetadot0_vec(k);

    if thetadot0 >= 0
        theta0 = gamma - alpha;
    else
        theta0 = gamma + alpha;
    end

    y0 = [theta0; thetadot0];
    double_support = 0;

    t0 = 0; tf = 25; dt = 0.01;
    tol = 1e-3;
    impact_counter = 0;

    while t0 < tf
        options = odeset('Events', @(t, y) impact_event(t, y, alpha, gamma), 'MaxStep', dt);
        [t, y, te, ye, ~] = ode45(@(t, y) dynamics(t, y, g, l, double_support), [t0 tf], y0, options);

        if ~isempty(te)
            impact_counter = impact_counter + 1;
            [y0, double_support] = impact_map(ye, alpha, g, l);
            t0 = te;
        else
            break;
        end
    end

    if double_support == 1
        colors(k, :) = [1 0 0]; % rosso = fermo
    else
        colors(k, :) = [0 0 1]; % blu = limit cycle
    end
end

% === Plot 2D Basin of Attraction ===
idx_red = colors(:,1) == 1;
idx_blue = colors(:,3) == 1;  
figure; hold on; grid on;
scatter(thetadot0_vec(idx_red), zeros(sum(idx_red),1), 10, 'r', 'filled');
scatter(thetadot0_vec(idx_blue), zeros(sum(idx_blue),1), 10, 'b', 'filled');
xlabel('$\dot{\theta}_0$ (rad/s)', 'Interpreter', 'latex');
ylabel('State classification', 'Interpreter', 'latex');
yticks([]);
title('Basin of Attraction of Rimless Wheel', 'Interpreter', 'latex');
legend({'Red: equilibrium', 'Blue: limit cycle'}, 'Location', 'best', 'Interpreter', 'latex');
set(gcf, 'PaperPositionMode', 'auto');   
print(gcf, 'eq_classification1.pdf', '-dpdf', '-bestfit');



%% 4b
clc; close all; clear all;

% Parametri
g = 9.81;
l_vec = [0.2, 0.5, 1, 2, 3, 5];
colors = rand(length(l_vec), 3);
alpha = pi/8;
gamma = 0.08;

thetadot0 = 0.95;

figure(1); hold on; grid on;
title('Phase Portraits for different l', 'Interpreter', 'latex');
xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');


figure(2)
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('State', 'Interpreter', 'latex');
title('Rimless Wheel Dynamics over Time', 'Interpreter', 'latex');
legend show;
set(legend, 'Interpreter', 'latex');


figure(3); clf; hold on; grid on;
title('Time intervals between impacts for different l', 'Interpreter', 'latex');
xlabel('Impact number', 'Interpreter', 'latex');
ylabel('Time interval $\Delta t$ (s)', 'Interpreter', 'latex');


for k = 1:length(l_vec)

    l = l_vec(k);
    if thetadot0 >= 0
        theta0 = gamma - alpha;
    else
        theta0 = gamma + alpha;
    end

    y0 = [theta0; thetadot0];
    double_support = 0;

    T = [];
    Y = [];
    impactTimes = []; 
    t0 = 0; tf = 25; dt = 0.01;

    while t0 < tf
        options = odeset('Events', @(t, y) impact_event(t, y, alpha, gamma), 'MaxStep', dt);
        [t, y, te, ye, ie] = ode45(@(t, y) dynamics(t, y, g, l, double_support), [t0 tf], y0, options);

        T = [T; t];
        Y = [Y; y];

        if ~isempty(te)
            impactTimes = [impactTimes; te(end)]; 

            [y0, double_support] = impact_map(ye(end,:)', alpha, g, l); % ye(end,:) trasposto colonna
            t0 = te(end);   % usa l'ultimo impatto
        else
            break;
        end
    end

    % === PLOT ===
    figure(1)
    subplot(3, 2, k)
    plot(Y(:,1), Y(:,2), 'Color', colors(k,:), 'LineWidth', 1);
    hold on;
    plot(Y(1,1), Y(1,2), '*', 'Color', colors(k,:), 'LineWidth', 1.5);
    grid on;
    xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
    ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');
    title(['Phase Portrait for $l= $' num2str(l)], 'Interpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');   
    print(gcf, 'rimless_limit_l.pdf', '-dpdf', '-bestfit');
    
    figure(2)
    subplot(3, 2, k)
    plot(T, Y(:,1), '-', 'Color', colors(k,:), 'LineWidth', 0.8);
    hold on;
    plot(T, Y(:,2), '--', 'Color', colors(k,:), 'LineWidth', 0.8);
    grid on;
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('State', 'Interpreter', 'latex');
    title(['Dynamics for $l = $' num2str(l)], 'Interpreter', 'latex');
    legend({'$\theta$', '$\dot{\theta}$'}, 'Interpreter', 'latex', 'Location', 'best'); 
    set(gcf, 'PaperPositionMode', 'auto');   
    print(gcf, 'rimless_dyn_l.pdf', '-dpdf', '-bestfit');

    figure(3)
    if length(impactTimes) > 1
        delta_t = diff(impactTimes);
        plot(1:length(delta_t), delta_t, '-o', 'Color', colors(k,:), 'DisplayName', ['$l=$' num2str(l)], 'LineWidth', 1);
        legend('Interpreter', 'latex', 'Location', 'best');
        set(gcf, 'PaperPositionMode', 'auto');   
        print(gcf, 'rimless_impact_l.pdf', '-dpdf', '-bestfit');
    else
        plot(1, NaN, 'x', 'Color', colors(k,:), 'DisplayName', ['$l=$' num2str(l)]);
        legend('Interpreter', 'latex', 'Location', 'best');
        set(gcf, 'PaperPositionMode', 'auto');   
        print(gcf, 'rimless_impact_l.pdf', '-dpdf', '-bestfit');
    end

end

%% 
clc; close all; clear all;

% Parametri
g = 9.81;
l = 1.0;
alpha_vec = [pi/32, pi/16, pi/12, pi/8, pi/6, pi/4];
gamma = 0.08;
colors = rand(length(alpha_vec), 3);
thetadot0 = 0.95;

figure(1); hold on; grid on;
title('Phase Portraits for different $\alpha$', 'Interpreter', 'latex');
xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');


figure(2)
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('State', 'Interpreter', 'latex');
title('Rimless Wheel Dynamics over Time', 'Interpreter', 'latex');
legend show;
set(legend, 'Interpreter', 'latex');


figure(3); clf; hold on; grid on;
title('Time intervals between impacts for different $\alpha$', 'Interpreter', 'latex');
xlabel('Impact number', 'Interpreter', 'latex');
ylabel('Time interval $\Delta t$ (s)', 'Interpreter', 'latex');


for k = 1:length(alpha_vec)

    alpha = alpha_vec(k);
    if thetadot0 >= 0
        theta0 = gamma - alpha;
    else
        theta0 = gamma + alpha;
    end

    y0 = [theta0; thetadot0];
    double_support = 0;

    T = [];
    Y = [];
    impactTimes = []; 
    t0 = 0; tf = 25; dt = 0.01;

    while t0 < tf
        options = odeset('Events', @(t, y) impact_event(t, y, alpha, gamma), 'MaxStep', dt);
        [t, y, te, ye, ie] = ode45(@(t, y) dynamics(t, y, g, l, double_support), [t0 tf], y0, options);

        T = [T; t];
        Y = [Y; y];

        if ~isempty(te)
            impactTimes = [impactTimes; te(end)]; 

            [y0, double_support] = impact_map(ye(end,:)', alpha, g, l); % ye(end,:) trasposto colonna
            t0 = te(end);   % usa l'ultimo impatto
        else
            break;
        end
    end

    % === PLOT ===
    figure(1)
    subplot(3, 2, k)
    plot(Y(:,1), Y(:,2), 'Color', colors(k,:), 'LineWidth', 1);
    hold on;
    plot(Y(1,1), Y(1,2), '*', 'Color', colors(k,:), 'LineWidth', 1.5);
    grid on;
    xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
    ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');
    title(['Phase Portrait for $\alpha= $' num2str(alpha)], 'Interpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');   
    print(gcf, 'rimless_limit_a.pdf', '-dpdf', '-bestfit');
    
    figure(2)
    subplot(3, 2, k)
    plot(T, Y(:,1), '-', 'Color', colors(k,:), 'LineWidth', 0.8);
    hold on;
    plot(T, Y(:,2), '--', 'Color', colors(k,:), 'LineWidth', 0.8);
    grid on;
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('State', 'Interpreter', 'latex');
    title(['Dynamics for $\alpha = $' num2str(alpha)], 'Interpreter', 'latex');
    legend({'$\theta$', '$\dot{\theta}$'}, 'Interpreter', 'latex', 'Location', 'best'); 
    set(gcf, 'PaperPositionMode', 'auto');   
    print(gcf, 'rimless_dyn_a.pdf', '-dpdf', '-bestfit');

    figure(3)
    if length(impactTimes) > 1
        delta_t = diff(impactTimes);
        plot(1:length(delta_t), delta_t, '-o', 'Color', colors(k,:), 'DisplayName', ['$\alpha=$' num2str(alpha)], 'LineWidth', 1);
        legend('Interpreter', 'latex', 'Location', 'best');
        set(gcf, 'PaperPositionMode', 'auto');   
        print(gcf, 'rimless_impact_a.pdf', '-dpdf', '-bestfit');
    else
        plot(1, NaN, 'x', 'Color', colors(k,:), 'DisplayName', ['$\alpha=$' num2str(alpha)]);
        legend('Interpreter', 'latex', 'Location', 'best');
        set(gcf, 'PaperPositionMode', 'auto');   
        print(gcf, 'rimless_impact_a.pdf', '-dpdf', '-bestfit');
    end

end


%%
%% 
clc; close all; clear all;

% Parametri
g = 9.81;
l = 1.0;
alpha =pi/8;
gamma_vec = [0.05, 0.1, 0.2, 0.5];
colors = rand(length(gamma_vec), 3);
thetadot0 = 0.95;

figure(1); hold on; grid on;
title('Phase Portraits for different $\gamma$', 'Interpreter', 'latex');
xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');


figure(2)
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('State', 'Interpreter', 'latex');
title('Rimless Wheel Dynamics over Time', 'Interpreter', 'latex');
legend show;
set(legend, 'Interpreter', 'latex');


figure(3); clf; hold on; grid on;
title('Time intervals between impacts for different $\gamma$', 'Interpreter', 'latex');
xlabel('Impact number', 'Interpreter', 'latex');
ylabel('Time interval $\Delta t$ (s)', 'Interpreter', 'latex');


for k = 1:length(gamma_vec)

    gamma = gamma_vec(k);
    if thetadot0 >= 0
        theta0 = gamma - alpha;
    else
        theta0 = gamma + alpha;
    end

    y0 = [theta0; thetadot0];
    double_support = 0;

    T = [];
    Y = [];
    impactTimes = []; 
    t0 = 0; tf = 25; dt = 0.01;

    while t0 < tf
        options = odeset('Events', @(t, y) impact_event(t, y, alpha, gamma), 'MaxStep', dt);
        [t, y, te, ye, ie] = ode45(@(t, y) dynamics(t, y, g, l, double_support), [t0 tf], y0, options);

        T = [T; t];
        Y = [Y; y];

        if ~isempty(te)
            impactTimes = [impactTimes; te(end)]; 

            [y0, double_support] = impact_map(ye(end,:)', alpha, g, l); % ye(end,:) trasposto colonna
            t0 = te(end);   % usa l'ultimo impatto
        else
            break;
        end
    end

    % === PLOT ===
    figure(1)
    subplot(2, 2, k)
    plot(Y(:,1), Y(:,2), 'Color', colors(k,:), 'LineWidth', 1);
    hold on;
    plot(Y(1,1), Y(1,2), '*', 'Color', colors(k,:), 'LineWidth', 1.5);
    grid on;
    xlabel('$\theta$ (rad)', 'Interpreter', 'latex');
    ylabel('$\dot{\theta}$ (rad/s)', 'Interpreter', 'latex');
    title(['Phase Portrait for $\gamma= $' num2str(gamma)], 'Interpreter', 'latex');
    set(gcf, 'PaperPositionMode', 'auto');   
    print(gcf, 'rimless_limit_gamma.pdf', '-dpdf', '-bestfit');
    
    figure(2)
    subplot(2, 2, k)
    plot(T, Y(:,1), '-', 'Color', colors(k,:), 'LineWidth', 0.8);
    hold on;
    plot(T, Y(:,2), '--', 'Color', colors(k,:), 'LineWidth', 0.8);
    grid on;
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('State', 'Interpreter', 'latex');
    title(['Dynamics for $\gamma = $' num2str(gamma)], 'Interpreter', 'latex');
    legend({'$\theta$', '$\dot{\theta}$'}, 'Interpreter', 'latex', 'Location', 'best'); 
    set(gcf, 'PaperPositionMode', 'auto');   
    print(gcf, 'rimless_dyn_gamma.pdf', '-dpdf', '-bestfit');

    figure(3)
    if length(impactTimes) > 1
        delta_t = diff(impactTimes);
        plot(1:length(delta_t), delta_t, '-o', 'Color', colors(k,:), 'DisplayName', ['$\gamma=$' num2str(gamma)], 'LineWidth', 1);
        legend('Interpreter', 'latex', 'Location', 'best');
        set(gcf, 'PaperPositionMode', 'auto');   
        print(gcf, 'rimless_impact_a.pdf', '-dpdf', '-bestfit');
    else
        plot(1, NaN, 'x', 'Color', colors(k,:), 'DisplayName', ['$\gamma=$' num2str(gamma)]);
        legend('Interpreter', 'latex', 'Location', 'best');
        set(gcf, 'PaperPositionMode', 'auto');   
        print(gcf, 'rimless_impact_gamma.pdf', '-dpdf', '-bestfit');
    end

end






function dydt = dynamics(~, y, g, l, ds)
    theta = y(1);
    thetadot = y(2);
    if (~ds)
        dtheta = thetadot;
        dthetadot = (g/l) * sin(theta);
    else
        dtheta = 0;
        dthetadot = 0;
    end
    dydt = [dtheta; dthetadot];
end

function [value, isterminal, direction] = impact_event(~, y, alpha,gamma)
    
    value = [y(1)-alpha-gamma; y(1)-gamma+alpha];% Trigger when theta = gamma+alpha
                                     %Trigger when theta = gamma-alpha
    isterminal = [1;1];         % Stop the integration
    direction = [1;-1];          % Detect only when increasing
end

function [yplus,ds] = impact_map(y_minus, alpha,g,l)%minus: before impact time; plus: after impact time
    if (y_minus(2)>=0)
        theta_plus = y_minus(1)-2*alpha;
    else
        theta_plus = y_minus(1)+2*alpha;
    end
    thetadot_plus = cos(2*alpha) * y_minus(2);
    if (thetadot_plus < 0.01*sqrt(g/l) && thetadot_plus >-0.01*sqrt(g/l)) 
        thetadot_plus = 0;
        ds = 1;
    else
        ds = 0;
    end
    yplus = [theta_plus; thetadot_plus];
end

