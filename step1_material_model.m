% Step 1: Silicone material models (hyperelastic)
% Silicone stiffness changes as it stretches, so we use hyperelastic laws instead of Young's modulus.
%   Neo-Hookean:  W = C10 (I1 - 3)
%   Yeoh:         W = C10 (I1 - 3) + C20 (I1 - 3)^2 + C30 (I1 - 3)^3
% Uniaxial nominal stress (incompressible):  P = 2 (lam - lam^-2) dW/dI1,  with I1 = lam^2 + 2/lam
% Parameter values (kPa) from the literature table in Runge et al., arXiv:2209.03831:
%   Ecoflex 00-30 C10 = 12.7, Ecoflex 00-50 C10 = 25, Dragon Skin 10 C10 = 42.5,
%   Elastosil M4601 neo-Hookean C10 = 131.2; Yeoh C10 = 110, C20 = 20 (Mosadegh et al.)
% Note: shear modulus mu = 2*C10. The Step 2-5 models use mu = 25 kPa (Ecoflex 00-30).
clear; clc; close all;

lam = linspace(1, 3, 200);
I1  = lam.^2 + 2./lam;
P_nh   = @(C10) 2*(lam - lam.^-2) .* C10;
P_yeoh = @(C10, C20, C30) 2*(lam - lam.^-2) .* (C10 + 2*C20*(I1 - 3) + 3*C30*(I1 - 3).^2);

mats = {'Ecoflex 00-30', 12.7; 'Ecoflex 00-50', 25; 'Dragon Skin 10', 42.5; 'Elastosil M4601', 131.2};

fig = figure('Color', 'w', 'Position', [100 100 1100 420]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on; grid on;
for i = 1:size(mats, 1)
    plot(lam, P_nh(mats{i,2}), 'LineWidth', 2, 'DisplayName', sprintf('%s (C10 = %.1f kPa)', mats{i,1}, mats{i,2}));
end
xlabel('Stretch \lambda'); ylabel('Nominal stress (kPa)');
title('Neo-Hookean, uniaxial tension'); legend('Location', 'northwest');

nexttile; hold on; grid on;
plot(lam, P_nh(131.2), 'LineWidth', 2, 'DisplayName', 'Elastosil neo-Hookean (C10 = 131.2)');
plot(lam, P_nh(110),   'LineWidth', 2, 'DisplayName', 'Elastosil neo-Hookean (C10 = 110)');
plot(lam, P_yeoh(110, 20, 0), 'LineWidth', 2, 'DisplayName', 'Elastosil Yeoh (C10 = 110, C20 = 20)');
xlabel('Stretch \lambda'); ylabel('Nominal stress (kPa)');
title('Where neo-Hookean and Yeoh disagree'); legend('Location', 'northwest');
exportgraphics(fig, 'step1_material_model.png', 'Resolution', 150);

% How far apart are neo-Hookean (C10 = 110) and Yeoh at typical finger strains?
fprintf('Stretch | Yeoh vs neo-Hookean (C10=110) stress difference\n');
for l = [1.1 1.3 1.5 2.0 2.5]
    I = l^2 + 2/l;
    nh = 2*(l - l^-2)*110;  ye = 2*(l - l^-2)*(110 + 2*20*(I - 3));
    fprintf('  %.1f   | %+5.1f %%\n', l, 100*(ye - nh)/nh);
end
fprintf('\nShear modulus mu = 2*C10 used in later steps:\n');
for i = 1:size(mats, 1), fprintf('  %-16s mu = %5.1f kPa\n', mats{i,1}, 2*mats{i,2}); end
