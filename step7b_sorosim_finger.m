function step7b_sorosim_finger()
% Step 7b: The soft finger itself in SoRoSim, compared with the Step 2 beam model.
% Written as a function so it does not touch variables in the MATLAB workspace.
%
% SoRoSim finger (built through the GUI, saved in finger_sorosim.mat):
%   100 x 20 x 17 mm solid rectangular rod, E = 76 kPa, fixed base, no gravity,
%   bending about y (order 2) + elongation, one straight cable 8 mm below the centre line.
% Pressure -> cable tension: matched by bending moment. In Step 2 the chamber pressure creates
%   M = P*Ac*yc (Ac = 192 mm^2, yc = 9 mm from the strain-limiting bottom layer).
%   The cable gives M = T*0.008, so T = P*Ac*yc/0.008.
% The two models differ on purpose:
%   Step 2  - hyperelastic, hollow chamber, bends about an inextensible bottom layer.
%   SoRoSim - linear elastic, SOLID section, bends about its centre (no strain-limiting layer).
% Real measurements in the lab decide which assumptions matter.

here = fileparts(mfilename('fullpath'));
base = fullfile(getenv('APPDATA'), 'MathWorks', 'MATLAB Add-Ons', 'Toolboxes', 'SoRoSim', 'Differentiable_SoRoSim');
if isempty(which('SorosimLinkage')), addpath(genpath(base)); end
addpath(fullfile(here, 'sorosim_shims'));
outdir = fullfile(here, 'sorosim_runs'); if ~exist(outdir, 'dir'), mkdir(outdir); end
old = cd(outdir);  c = onCleanup(@() cd(old));

D = load(fullfile(here, 'finger_sorosim.mat'));  Fg = D.F;
fprintf('SoRoSim finger: %d DOF, %d actuator(s)\n', Fg.ndof, Fg.nact);

Ac = 0.016*0.012;  yc = 0.009;  d = 0.008;
P  = 0:1e3:15e3;                          % pressure sweep [Pa]
T  = P*Ac*yc/d;                           % equivalent cable tension [N]

thSoRo = zeros(size(P));  shapes = cell(size(P));
q = zeros(Fg.ndof, 1);
for k = 1:numel(P)
    if P(k) > 0
        evalc('[q,~,~] = Fg.statics(q, T(k));');
        close all;
    end
    g = Fg.FwdKinematics(q);  shapes{k} = g;
    gt = g(end-3:end, :);
    thSoRo(k) = rad2deg(atan2(gt(3,1), gt(1,1)));       % rotation about y (bend in x-z plane)
end
thSoRo = abs(thSoRo);

% Step 2 beam model (hyperelastic, inextensible bottom layer)
mu = 25.4e3; L = 0.10; b = 0.020; wall = 0.002; yb = 0.003; yt = 0.015; H = 0.017;
s    = @(lam) mu*(lam - 1./lam.^2);
bw   = @(y) b*(y < yb | y > yt) + 2*wall*(y >= yb & y <= yt);
Mres = @(kk) integral(@(y) s(1 + kk*y).*y.*bw(y), 0, H, 'Waypoints', [yb yt]);
thBeam = zeros(size(P));
for k = 2:numel(P)
    thBeam(k) = rad2deg(fzero(@(kk) P(k)*Ac*yc - Mres(kk), [0 500])*L);
end

fprintf('\nPressure | cable T | SoRoSim bend | Step 2 bend\n');
for k = 1:numel(P)
    fprintf('%5.1f kPa | %5.2f N | %8.1f deg  | %7.1f deg\n', P(k)/1e3, T(k), thSoRo(k), thBeam(k));
end

fig = figure('Color', 'w', 'Position', [100 100 1150 430]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on; grid on; axis equal;
for k = 1:3:numel(P)
    g = shapes{k};  n = size(g,1)/4;
    plot(g(4*(0:n-1)+1, 4)*1e3, g(4*(0:n-1)+3, 4)*1e3, 'LineWidth', 2, 'DisplayName', sprintf('%g kPa', P(k)/1e3));
end
xlabel('x (mm)'); ylabel('z (mm)'); title('SoRoSim finger shape (centre line)'); legend('Location', 'southwest');
nexttile; hold on; grid on;
plot(P/1e3, thSoRo, 'o-', 'LineWidth', 2, 'DisplayName', 'SoRoSim (solid, linear elastic)');
plot(P/1e3, thBeam, 's-', 'LineWidth', 2, 'DisplayName', 'Step 2 (hollow, hyperelastic)');
xlabel('Pressure (kPa)'); ylabel('Tip bend (deg)'); title('Two models, same finger'); legend('Location', 'northwest');
pngFile = fullfile(here, 'step7b_sorosim_finger.png');
try
    exportgraphics(fig, pngFile, 'Resolution', 150);
catch
    % The PNG is locked (open in a viewer, or OneDrive is syncing it) - save under a new name instead
    pngFile = fullfile(here, ['step7b_sorosim_finger_' char(datetime('now','Format','HHmmss')) '.png']);
    exportgraphics(fig, pngFile, 'Resolution', 150);
end
fprintf('\nPlot saved: %s\n', pngFile);
end
