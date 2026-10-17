% Step 3: Parameter sweep - how geometry and material change the bend angle
% Same first-order beam model as Step 2. Values are rough starting estimates.
clear; clc; close all;

base.mu = 25e3; base.L = 0.10; base.b = 0.020; base.H = 0.017;
base.yb = 0.003; base.yt = 0.015; base.wall = 0.002;
P = linspace(0, 25e3, 40);

fig = figure('Color', 'w', 'Position', [100 100 1200 380]);
tiledlayout(1, 3, 'TileSpacing', 'compact');

% 1) Side wall thickness
nexttile; hold on; grid on;
for w = [1 2 3 4]*1e-3
    p = base; p.wall = w;
    plot(P/1e3, bendCurve(p, P), 'LineWidth', 2, 'DisplayName', sprintf('%g mm', w*1e3));
end
title('Side wall thickness'); xlabel('Pressure (kPa)'); ylabel('Bend angle (deg)');
legend('Location', 'northwest'); ylim([0 360]); yline(180, '--', 'HandleVisibility', 'off');

% 2) Bottom (strain-limiting side) thickness
nexttile; hold on; grid on;
for yb = [1 2 3 5]*1e-3
    p = base; p.yb = yb;
    plot(P/1e3, bendCurve(p, P), 'LineWidth', 2, 'DisplayName', sprintf('%g mm', yb*1e3));
end
title('Bottom layer thickness'); xlabel('Pressure (kPa)');
legend('Location', 'northwest'); ylim([0 360]); yline(180, '--', 'HandleVisibility', 'off');

% 3) Material (shear modulus mu = 2*C10, literature values - see step1_material_model.m)
nexttile; hold on; grid on;
mats = {'Ecoflex 00-30', 25.4e3; 'Ecoflex 00-50', 50e3; 'Dragon Skin 10', 85e3; 'Elastosil M4601', 262.4e3};
for i = 1:size(mats, 1)
    p = base; p.mu = mats{i, 2};
    plot(P/1e3, bendCurve(p, P), 'LineWidth', 2, 'DisplayName', mats{i, 1});
end
title('Silicone material'); xlabel('Pressure (kPa)');
legend('Location', 'northwest'); ylim([0 360]); yline(180, '--', 'HandleVisibility', 'off');

exportgraphics(fig, 'step3_param_sweep.png', 'Resolution', 150);

% Pressure needed for a 90 deg bend, per material
fprintf('Pressure for 90 deg bend:\n');
for i = 1:size(mats, 1)
    p = base; p.mu = mats{i, 2};
    th = bendCurve(p, linspace(0, 150e3, 150));
    Pq = linspace(0, 150e3, 150);
    idx = find(th >= 90, 1);
    if isempty(idx), fprintf('  %-15s > 150 kPa\n', mats{i,1});
    else, fprintf('  %-15s %.1f kPa\n', mats{i,1}, Pq(idx)/1e3); end
end

function theta = bendCurve(p, P)
    Ac = (p.b - 2*p.wall) * (p.yt - p.yb);
    yc = (p.yb + p.yt) / 2;
    s    = @(lam) p.mu * (lam - 1 ./ lam.^2);
    bw   = @(y) p.b * (y < p.yb | y > p.yt) + 2*p.wall * (y >= p.yb & y <= p.yt);
    Mres = @(k) integral(@(y) s(1 + k*y) .* y .* bw(y), 0, p.H, 'Waypoints', [p.yb p.yt]);
    theta = zeros(size(P));
    for i = 2:numel(P)
        f = @(k) P(i)*Ac*yc - Mres(k);
        if f(500) > 0, theta(i) = NaN; continue; end   % beyond model range
        k = fzero(f, [0 500]);
        theta(i) = rad2deg(k * p.L);
    end
end
