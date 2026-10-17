% Step 4: Pressure dynamics - pump/valve -> air chamber -> finger bend over time
% Air is treated as an ideal gas at constant temperature (isothermal).
% Valve = small orifice (standard compressible-flow orifice equation).
% Chamber volume grows with pressure (compliance C).
% Finger motion = static bend from Step 2 + a springy mass-damper lag.
% All parameter values are rough estimates to be calibrated in the lab.
clear; clc; close all;

% ---- Air and supply ----
R = 287; T = 293; k = 1.4; patm = 101325;     % [J/kgK], [K], [-], [Pa]
Psup_g = 10e3;                                % supply pressure, gauge [Pa]
Psup = patm + Psup_g;
Cd = 0.6;                                     % valve discharge coefficient

% ---- Chamber (same geometry as Step 2) ----
L = 0.10; b = 0.020; wall = 0.002; yb = 0.003; yt = 0.015;
V0 = (b - 2*wall) * (yt - yb) * L;            % ~19 mL at rest
C  = V0 / 25e3;                               % volume roughly doubles by 25 kPa (guess)

% ---- Finger mechanics (guesses) ----
wn = 2*pi*4;  zeta = 0.15;                    % 4 Hz natural frequency, light damping

% Static pressure -> bend map from the Step 2 model
p.mu = 25e3; p.L = L; p.b = b; p.H = 0.017; p.yb = yb; p.yt = yt; p.wall = wall;
Ptab  = linspace(0, Psup_g*1.05, 40);
thtab = bendCurve(p, Ptab);

diams = [0.3 0.6 1.2] * 1e-3;                 % valve orifice diameters [m]
phases = [0 1; 1 2; 2 3];                     % inflate, hold, vent [s]
valve  = [1 0; 0 0; 0 1];                     % [inlet exhaust] open flags

fig = figure('Color', 'w', 'Position', [100 100 1100 420]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on');

fprintf('Orifice  | 90%% pressure rise | peak bend | final bend (inflate)\n');
for d = diams
    A = pi * (d/2)^2;
    x = [patm; 0; 0];  tAll = []; xAll = [];
    for ph = 1:3
        f = @(t, x) rhs(x, valve(ph,:), A, Cd, Psup, patm, R, T, k, V0, C, Ptab, thtab, wn, zeta);
        opts = odeset('RelTol', 1e-6, 'AbsTol', [1 1e-6 1e-6]);
        [t, X] = ode15s(f, phases(ph,:), x, opts);
        tAll = [tAll; t]; xAll = [xAll; X]; x = X(end, :)'; %#ok<AGROW>
    end
    Pg = (xAll(:,1) - patm) / 1e3;
    lbl = sprintf('%.1f mm valve', d*1e3);
    plot(ax1, tAll, Pg, 'LineWidth', 2, 'DisplayName', lbl);
    plot(ax2, tAll, xAll(:,2), 'LineWidth', 2, 'DisplayName', lbl);

    inflate = tAll <= 1;
    i90 = find(Pg(inflate) >= 0.9*Psup_g/1e3, 1);
    if isempty(i90), t90 = NaN; else, t90 = tAll(i90); end   % NaN = never reached 90% in 1 s
    fprintf('%.1f mm   | %8.3f s        | %6.1f deg | %6.1f deg\n', ...
        d*1e3, t90, max(xAll(inflate,2)), xAll(find(inflate,1,'last'),2));
end

for ax = [ax1 ax2], xline(ax, [1 2], ':', {'hold', 'vent'}, 'HandleVisibility', 'off'); end
xlabel(ax1, 'Time (s)'); ylabel(ax1, 'Chamber pressure (kPa, gauge)');
title(ax1, 'Pressure: inflate 0-1 s, hold, vent 2-3 s'); legend(ax1, 'Location', 'east');
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'Bend angle (deg)');
title(ax2, 'Finger bend over time'); legend(ax2, 'Location', 'east');
exportgraphics(fig, 'step4_pressure_dynamics.png', 'Resolution', 150);

% ------------------------------------------------------------------------
function dx = rhs(x, open, A, Cd, Psup, patm, R, T, k, V0, C, Ptab, thtab, wn, zeta)
    pc = x(1); th = x(2); om = x(3);
    mdot = open(1)*orifice(Psup, pc, A, Cd, R, T, k) - open(2)*orifice(pc, patm, A, Cd, R, T, k);
    V  = max(V0 + C*(pc - patm), V0);
    dp = R*T*mdot / (V + pc*C);                       % isothermal chamber with growing volume
    thss = interp1(Ptab, thtab, min(max(pc - patm, 0), Ptab(end)));
    dom = wn^2*(thss - th) - 2*zeta*wn*om;            % springy finger lags the pressure
    dx = [dp; om; dom];
end

function m = orifice(pu, pd, A, Cd, R, T, k)
    % Compressible flow through an orifice; positive from pu to pd
    s = 1;
    if pd > pu, [pu, pd] = deal(pd, pu); s = -1; end
    r = pd / pu; rc = (2/(k+1))^(k/(k-1));
    if r <= rc   % choked
        m = Cd*A*pu*sqrt(k/(R*T)) * (2/(k+1))^((k+1)/(2*(k-1)));
    else
        m = Cd*A*pu*sqrt(2*k/(R*T*(k-1)) * (r^(2/k) - r^((k+1)/k)));
    end
    m = s*m;
end

function theta = bendCurve(p, P)
    Ac = (p.b - 2*p.wall) * (p.yt - p.yb);  yc = (p.yb + p.yt)/2;
    s    = @(lam) p.mu * (lam - 1 ./ lam.^2);
    bw   = @(y) p.b*(y < p.yb | y > p.yt) + 2*p.wall*(y >= p.yb & y <= p.yt);
    Mres = @(k) integral(@(y) s(1 + k*y) .* y .* bw(y), 0, p.H, 'Waypoints', [p.yb p.yt]);
    theta = zeros(size(P));
    for i = 2:numel(P)
        k = fzero(@(k) P(i)*Ac*yc - Mres(k), [0 500]);
        theta(i) = rad2deg(k * p.L);
    end
end
