% Step 5: Closed-loop control of the soft finger's bend angle
% Compares two controllers on a "real" finger that is 20% stiffer than our model thinks:
%   A) Pressure-only control: look up the pressure our model says gives the target angle,
%      then hold that pressure (no angle sensor).
%   B) Cascade angle PID: measure the bend angle (flex sensor / IMU). An outer PID on angle
%      sets the target pressure (starting from the model's guess), and a fast inner loop
%      drives the valves to reach that pressure. (A single PID straight from angle to valves
%      oscillates badly here: air pressure integrates valve flow and the finger is springy.)
% Valves: one inlet + one exhaust solenoid (1.2 mm), PWM duty treated as an averaged opening.
% Controller runs at 200 Hz with sensor noise, like an ESP32 loop would.
clear; clc; close all;

P.R = 287; P.T = 293; P.k = 1.4; P.patm = 101325;
P.Psup = P.patm + 20e3;  P.Cd = 0.6;  P.A = pi*(1.2e-3/2)^2;
L = 0.10; b = 0.020; wall = 0.002; yb = 0.003; yt = 0.015; H = 0.017;
P.V0 = (b - 2*wall)*(yt - yb)*L;  P.C = P.V0/25e3;
P.wn = 2*pi*4;  P.zeta = 0.15;

model = struct('mu', 25e3, 'L', L, 'b', b, 'H', H, 'yb', yb, 'yt', yt, 'wall', wall);
real  = model;  real.mu = 30e3;                     % the actual finger is 20% stiffer
P.Ptab = linspace(0, 20e3, 60);
P.thModel = bendCurve(model, P.Ptab);
P.thReal  = bendCurve(real,  P.Ptab);

dt = 2e-4; t = 0:dt:4.5;
ref = 60*(t >= 0.2) + 40*(t >= 1.7) - 70*(t >= 3.2);   % 0 -> 60 -> 100 -> 30 deg

% ---- Tune each controller with a small grid search ----
bestA = []; JA = inf;
for kpi = [0.1 0.2 0.4 0.8]          % inner pressure loop: valve command per kPa error
    [th, ~, ~] = simulate(1, [kpi 0 0 0 0.02], t, ref, P);
    J = cost(t, th, ref);
    if J < JA, JA = J; bestA = [kpi 0 0 0 0.02]; end
end
bestB = []; JB = inf;
for kp = [0.02 0.05 0.1]             % outer angle loop: kPa per deg error
    for ki = [0.1 0.3 1]
        for kd = [0 0.002 0.005]
            g = [bestA(1) kp ki kd 0.02];
            [th, ~, ~] = simulate(2, g, t, ref, P);
            J = cost(t, th, ref);
            if J < JB, JB = J; bestB = g; end
        end
    end
end

[thA, pA, uA] = simulate(1, bestA, t, ref, P);
[thB, pB, uB] = simulate(2, bestB, t, ref, P);

fprintf('Inner pressure loop: Kp=%.2f per kPa\n', bestA(1));
fprintf('Outer angle PID:     Kp=%.3f kPa/deg, Ki=%.3f kPa/(deg s), Kd=%.4f kPa s/deg\n\n', bestB(2), bestB(3), bestB(4));
seg = [0.2 1.7 60; 1.7 3.2 100; 3.2 4.5 30];
fprintf('Step        | Controller     | Overshoot | Settle (+-2 deg) | Final error\n');
names = {'Pressure-only', 'Angle PID'};  TH = {thA, thB};
for s = 1:3
    for c = 1:2
        [os, ts, fe] = stepMetrics(t, TH{c}, seg(s,:));
        fprintf('to %3d deg  | %-14s | %6.1f deg | %10.2f s     | %+6.1f deg\n', ...
            seg(s,3), names{c}, os, ts, fe);
    end
end

fig = figure('Color', 'w', 'Position', [100 100 1150 720]);
tiledlayout(3, 1, 'TileSpacing', 'compact');
nexttile; hold on; grid on;
plot(t, ref, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Target');
plot(t, thA, 'LineWidth', 2, 'DisplayName', 'Pressure-only (no angle sensor)');
plot(t, thB, 'LineWidth', 2, 'DisplayName', 'Angle PID (with bend sensor)');
ylabel('Bend (deg)'); title('Bend angle tracking on a finger 20% stiffer than the model');
legend('Location', 'northeast');
nexttile; hold on; grid on;
plot(t, pA/1e3, 'LineWidth', 1.5, 'DisplayName', 'Pressure-only');
plot(t, pB/1e3, 'LineWidth', 1.5, 'DisplayName', 'Angle PID');
ylabel('Pressure (kPa)'); legend('Location', 'northeast');
nexttile; hold on; grid on;
plot(t, uB, 'LineWidth', 1.2);
ylabel('Valve command'); xlabel('Time (s)');
title('Angle PID valve command (+ = inlet PWM duty, - = exhaust PWM duty)');
exportgraphics(fig, 'step5_pid_control.png', 'Resolution', 150);

% ======================================================================
function [th, pg, u] = simulate(mode, g, t, ref, P)
    rng(1);                                  % same sensor noise for every run
    dt = t(2) - t(1);  N = numel(t);  Ts = 0.005;
    x = [P.patm; 0; 0];  th = zeros(1, N);  pg = th;  u = th;
    uc = 0;  I = 0;  Dm = 0;  prev = 0;  nextT = 0;  e = 0;
    Fmodel = griddedInterpolant(P.thModel, P.Ptab);   % angle -> pressure (what we believe)
    Freal  = griddedInterpolant(P.Ptab, P.thReal);    % pressure -> angle (what really happens)
    for i = 1:N
        if t(i) >= nextT - 1e-12                       % controller tick (200 Hz)
            nextT = nextT + Ts;
            thm = x(2) + 0.5*randn;                    % angle sensor, 0.5 deg noise
            pm  = (x(1) - P.patm)/1e3 + 0.05*randn;    % pressure sensor, kPa
            pff = Fmodel(min(ref(i), P.thModel(end)))/1e3;   % model's pressure guess [kPa]
            if mode == 1
                pref = pff;                                    % trust the model, no angle sensor
            else
                e = ref(i) - thm;
                Dm = Dm + (Ts/(g(5) + Ts))*((thm - prev)/Ts - Dm);   % filtered derivative
                prev = thm;
                Inew = I + e*Ts;
                pref = pff + g(2)*e + g(3)*Inew - g(4)*Dm;
                if pref > 0 && pref < 19, I = Inew; end            % anti-windup
                pref = min(max(pref, 0), 19);
            end
            uc = max(min(g(1)*(pref - pm), 1), -1);            % inner pressure loop -> valves
        end
        pc = x(1);
        mdot = orifice(P.Psup, pc, max(uc,0)*P.A, P) - orifice(pc, P.patm, max(-uc,0)*P.A, P);
        V  = max(P.V0 + P.C*(pc - P.patm), P.V0);
        dp = P.R*P.T*mdot/(V + pc*P.C);
        thss = Freal(min(max(pc - P.patm, 0), P.Ptab(end)));
        dom = P.wn^2*(thss - x(2)) - 2*P.zeta*P.wn*x(3);
        if x(2) < -5, dom = dom + 20*P.wn^2*(-5 - x(2)); end   % finger can't bend far backwards
        x = x + dt*[dp; x(3); dom];
        th(i) = x(2);  pg(i) = x(1) - P.patm;  u(i) = uc;
    end
end

function J = cost(t, th, ref)
    dt = t(2) - t(1);
    seg = [0.2 1.7 60; 1.7 3.2 100; 3.2 4.5 30];
    os = 0;
    for s = 1:3, os = os + stepMetrics(t, th, seg(s,:)); end
    J = sum(abs(ref - th))*dt + 0.5*os;
end

function [os, ts, fe] = stepMetrics(t, th, seg)
    idx = t >= seg(1) & t < seg(2);  tt = t(idx);  y = th(idx);  target = seg(3);
    start = th(find(idx, 1) - 1);
    if target >= start, os = max(max(y) - target, 0); else, os = max(target - min(y), 0); end
    out = find(abs(y - target) > 2, 1, 'last');
    if isempty(out), ts = 0; elseif out == numel(y), ts = NaN; else, ts = tt(out + 1) - seg(1); end
    fe = mean(y(tt > seg(2) - 0.2)) - target;
end

function m = orifice(pu, pd, A, P)
    if A <= 0, m = 0; return; end
    s = 1;  if pd > pu, [pu, pd] = deal(pd, pu); s = -1; end
    k = P.k;  r = pd/pu;  rc = (2/(k+1))^(k/(k-1));
    if r <= rc
        m = P.Cd*A*pu*sqrt(k/(P.R*P.T))*(2/(k+1))^((k+1)/(2*(k-1)));
    else
        m = P.Cd*A*pu*sqrt(2*k/(P.R*P.T*(k-1))*(r^(2/k) - r^((k+1)/k)));
    end
    m = s*m;
end

function theta = bendCurve(p, P)
    Ac = (p.b - 2*p.wall)*(p.yt - p.yb);  yc = (p.yb + p.yt)/2;
    s    = @(lam) p.mu*(lam - 1./lam.^2);
    bw   = @(y) p.b*(y < p.yb | y > p.yt) + 2*p.wall*(y >= p.yb & y <= p.yt);
    Mres = @(k) integral(@(y) s(1 + k*y).*y.*bw(y), 0, p.H, 'Waypoints', [p.yb p.yt]);
    theta = zeros(size(P));
    for i = 2:numel(P)
        f = @(k) P(i)*Ac*yc - Mres(k);
        if f(500) > 0, theta(i) = NaN; continue; end
        theta(i) = rad2deg(fzero(f, [0 500])*p.L);
    end
    last = find(~isnan(theta), 1, 'last');  theta(last+1:end) = theta(last);
end
