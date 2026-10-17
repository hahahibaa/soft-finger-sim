% Step 6: Grip force - how hard does the finger push once it touches an object?
% The finger bends freely until it reaches the contact angle theta_c, then the object blocks it.
% Extra pressure then turns into force at the fingertip instead of more bending.
% Lumped estimate (moment balance at the finger base, constant-curvature shape at contact):
%   F_tip = ( P*Ac*yc  -  M_rubber(curvature at contact) ) / L
% This is a first-order estimate: real fingers change shape under load, so expect the true
% blocked force to differ - measuring it with a load cell in the lab is the calibration step.
clear; clc; close all;

L = 0.10; b = 0.020; wall = 0.002; yb = 0.003; yt = 0.015; H = 0.017;
Ac = (b - 2*wall)*(yt - yb);  yc = (yb + yt)/2;
bw = @(y) b*(y < yb | y > yt) + 2*wall*(y >= yb & y <= yt);
Mres = @(mu, k) integral(@(y) mu*((1 + k*y) - (1 + k*y).^-2).*y.*bw(y), 0, H, 'Waypoints', [yb yt]);
Ftip = @(mu, P, thc) max(0, (P*Ac*yc - Mres(mu, deg2rad(thc)/L)) / L);

P = linspace(0, 100e3, 201);
fig = figure('Color', 'w', 'Position', [100 100 1100 420]);
tiledlayout(1, 2, 'TileSpacing', 'compact');

mats = {'Ecoflex 00-30', 25.4e3; 'Dragon Skin 10', 85e3};
for m = 1:2
    nexttile; hold on; grid on;
    for thc = [0 30 60 90]
        F = arrayfun(@(p) Ftip(mats{m,2}, p, thc), P);
        plot(P/1e3, F, 'LineWidth', 2, 'DisplayName', sprintf('contact at %d deg', thc));
    end
    xlabel('Pressure (kPa)'); ylabel('Fingertip force (N)');
    title(sprintf('%s finger', mats{m,1})); legend('Location', 'northwest');
end
exportgraphics(fig, 'step6_grip_force.png', 'Resolution', 150);

% ---- Example: 3-finger gripper holding an egg by friction ----
m_egg = 0.060; g = 9.81; mu_f = 0.5; n = 3; SF = 2;        % 60 g egg, silicone-shell friction ~0.5, safety factor 2
Fneed = SF*m_egg*g/(mu_f*n);                               % normal force needed per finger
fprintf('Egg (60 g), 3 fingers, friction 0.5, safety x2: each finger must push %.2f N\n\n', Fneed);
fprintf('Material        | contact | pressure to hold egg | force at 20 kPa\n');
for m = 1:2
    for thc = [30 60]
        F = arrayfun(@(p) Ftip(mats{m,2}, p, thc), P);
        i = find(F >= Fneed, 1);
        if isempty(i), ptxt = ' > 100 kPa'; else, ptxt = sprintf('%7.1f kPa', P(i)/1e3); end
        fprintf('%-15s | %3d deg | %s          | %.2f N\n', mats{m,1}, thc, ptxt, Ftip(mats{m,2}, 20e3, thc));
    end
end
fprintf('\nFor scale: a raw hen egg typically needs tens of newtons to crack, so these forces are gentle.\n');
fprintf('Pure friction pinch needs a lot of pressure; a wrapping (caging) grasp that cups the object\n');
fprintf('below its widest point needs far less force, which is how soft grippers usually hold eggs.\n');
