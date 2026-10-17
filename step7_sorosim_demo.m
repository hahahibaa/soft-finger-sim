% Step 7: First run of SoRoSim (Cosserat-rod soft robot simulator) on this laptop
% Uses the toolbox's built-in "SoftManipulator" example: one tapered soft rod (0.9 m, E = 1 MPa)
% with 6 internal cables. We pull cable 1 with increasing tension and record how the rod bends.
%
% Why this matters for the pneumatic finger: a pressurised chamber bends the finger with an
% internal moment M = P*Ac*yc. A cable at offset yc with tension T = P*Ac gives the same moment,
% so cable-actuated SoRoSim models are the quickest way to approximate a pneumatic finger
% (the axial force differs: pressure stretches the finger, a cable compresses it).
%
% Setup note: SoRoSim's prebuilt MEX files were compiled with a different MATLAB Coder version,
% so the folder 'sorosim_shims' redirects them to the toolbox's own MATLAB code (slower, same maths).
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
base = fullfile(getenv('APPDATA'), 'MathWorks', 'MATLAB Add-Ons', 'Toolboxes', 'SoRoSim', 'Differentiable_SoRoSim');
if isempty(which('SorosimLinkage')), addpath(genpath(base)); end
addpath(fullfile(here, 'sorosim_shims'));                   % must be first on the path
outdir = fullfile(here, 'sorosim_runs'); if ~exist(outdir, 'dir'), mkdir(outdir); end
old = cd(outdir);  cleanup = onCleanup(@() cd(old));       % SoRoSim saves .mat files to the current folder

ex = dir(fullfile(base, '**', 'SoftManipulator.mat'));
S  = load(fullfile(ex(1).folder, ex(1).name));  Rod = S.S1;

T = [0 0.25 0.5 1 2 3 4];             % cable-1 tension [N]
tipAngle = zeros(size(T));  shapes = cell(size(T));
q = zeros(Rod.ndof, 1);
for k = 1:numel(T)
    u = zeros(Rod.nact, 1);  u(1) = T(k);
    if T(k) > 0
        evalc('[q,~,~] = Rod.statics(q, u);');             % previous solution as the initial guess
        close all;
    end
    g = Rod.FwdKinematics(q);
    pts = g(4:4:end, :);  shapes{k} = g;                    % every 4th row block holds a frame
    gt = g(end-3:end, :);
    tipAngle(k) = rad2deg(atan2(-gt(2,1), gt(1,1)));
    fprintf('T = %4.2f N  ->  tip rotation %6.2f deg, tip at (x, y) = (%.3f, %.3f) m\n', T(k), tipAngle(k), gt(1,4), gt(2,4));
end

fig = figure('Color', 'w', 'Position', [100 100 1100 420]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on; grid on; axis equal;
for k = 1:numel(T)
    g = shapes{k};  n = size(g,1)/4;
    x = g(4*(0:n-1)+1, 4);  y = g(4*(0:n-1)+2, 4);
    plot(x, y, 'LineWidth', 2, 'DisplayName', sprintf('T = %.2f N', T(k)));
end
xlabel('x (m)'); ylabel('y (m)'); title('SoRoSim: rod shape as cable 1 is pulled'); legend('Location', 'southwest');
nexttile; hold on; grid on;
plot(T, tipAngle, 'o-', 'LineWidth', 2);
lin = T*(tipAngle(2)/T(2));
plot(T, lin, '--', 'LineWidth', 1.2);
legend('SoRoSim (geometrically exact rod)', 'Straight line through the smallest load', 'Location', 'northwest');
xlabel('Cable tension (N)'); ylabel('Tip rotation (deg)'); title('Bending vs tension');
exportgraphics(fig, fullfile(here, 'step7_sorosim_demo.png'), 'Resolution', 150);
