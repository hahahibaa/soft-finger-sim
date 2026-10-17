% Step 2: Pneumatic finger - pressure -> bend angle (first-order beam model)
% Strain-limited bottom layer at y = 0, chamber between yb and yt.
% Material: neo-Hookean silicone. Values are rough starting estimates.
clear; clc; close all;

mu   = 25e3;              % Ecoflex 00-30 shear modulus [Pa] (approx)
L    = 0.10;              % finger length [m]
b    = 0.020;             % width [m]
H    = 0.017;             % total height [m]
yb   = 0.003;             % chamber bottom height [m]
yt   = 0.015;             % chamber top height [m]
wall = 0.002;             % side wall thickness [m]

Ac = (b - 2*wall) * (yt - yb);   % chamber cross-section area [m^2]
yc = (yb + yt) / 2;              % lever arm of pressure force [m]

s    = @(lam) mu * (lam - 1 ./ lam.^2);                         % neo-Hookean uniaxial stress
bw   = @(y) b * (y < yb | y > yt) + 2*wall * (y >= yb & y <= yt); % material width at height y
Mres = @(k) integral(@(y) s(1 + k*y) .* y .* bw(y), 0, H, 'Waypoints', [yb yt]);

P     = linspace(0, 25e3, 50);   % pressure sweep [Pa]
theta = zeros(size(P));
for i = 2:numel(P)
    k = fzero(@(k) P(i)*Ac*yc - Mres(k), [0 500]);   % curvature [1/m]
    theta(i) = rad2deg(k * L);
end

fig = figure('Color', 'w');
plot(P/1e3, theta, 'LineWidth', 2); grid on;
xlabel('Pressure (kPa)'); ylabel('Bend angle (deg)');
title('Soft pneumatic finger: bend angle vs pressure');
exportgraphics(fig, 'step2_bend_vs_pressure.png', 'Resolution', 150);

fprintf('Bend at 10 kPa: %.1f deg\n', interp1(P, theta, 10e3));
fprintf('Bend at 20 kPa: %.1f deg\n', interp1(P, theta, 20e3));
fprintf('Bend at 25 kPa: %.1f deg\n', theta(end));
