# Soft pneumatic finger – simulation notes

Simulation-first groundwork for the December internship (soft robotic gripper lab).
All numbers are from simple models with literature or guessed parameters. They show **trends and design trade-offs**, not final values. The lab work is to measure the real finger and calibrate these models.

## Baseline finger (used in every step unless stated)
- Length 100 mm, width 20 mm, height 17 mm
- One air chamber 12 mm tall, 2 mm side walls, 3 mm strain-limiting bottom layer
- Silicone: Ecoflex 00-30, shear modulus about 25 kPa

## Files and results

| Step | Script | What it shows |
|---|---|---|
| 1 | `step1_material_model.m` | Silicone stress-stretch curves. Neo-Hookean is fine up to about 30% stretch (8% error vs Yeoh), but off by 21% at 50% stretch and 73% at 100%. Literature C10 values: Ecoflex 00-30 12.7 kPa, 00-50 25 kPa, Dragon Skin 10 42.5 kPa, Elastosil 131 kPa. |
| 2 | `step2_pressure_bend.m` | Pressure → bend angle. About 111° at 10 kPa. The model over-predicts at large bends (more than 180° past about 16 kPa), because it ignores wall ballooning. |
| 3 | `step3_param_sweep.m` | Pressure for a 90° bend: Ecoflex 00-30 9 kPa, 00-50 17 kPa, Dragon Skin 10 29 kPa, Elastosil 89 kPa. Thinner side walls bend more; bottom-layer thickness barely matters in this model (needs checking in the lab). |
| 4 | `step4_pressure_dynamics.m` | Valve size vs speed. A 1.2 mm valve reaches 90% pressure in 0.16 s but overshoots to 130° and wobbles; 0.6 mm takes 0.64 s and is smooth; 0.3 mm is too slow. |
| 5 | `step5_pid_control.m` | Closed loop on a finger 20% stiffer than the model. Pressure-only control ends up to 20° off target. Two-loop PID with a bend sensor (angle → target pressure → valves) is within 1.5° with almost no overshoot. A single PID from angle straight to valves oscillated badly. |
| 6 | `step6_grip_force.m` | Fingertip force once the finger touches an object. Only about 0.2–0.3 N at 20 kPa. Pinching a 60 g egg by friction (3 fingers) needs about 50 kPa (Ecoflex) or 57–66 kPa (Dragon Skin 10). A wrapping/cupping grasp needs far less force. |
| 7 | `step7_sorosim_demo.m` | SoRoSim works on this laptop. A cable-driven soft rod bends 14.5° per newton of tension, perfectly linearly (SoRoSim uses linear elasticity, unlike the hyperelastic model in Step 2). |
| 7b | `step7b_sorosim_finger.m` | The actual finger (100 × 20 × 17 mm, E = 76 kPa) in SoRoSim, built in the GUI and saved as `finger_sorosim.mat`. Pressure is mapped to cable tension by equal bending moment (10 kPa ≈ 2.2 N). SoRoSim bends 11.5° per kPa, linearly. Step 2 starts about 25% lower (8.9° at 1 kPa) and catches up by about 12 kPa (136–138° for both). Two very different models agreeing this closely is a good sign; the lab data decides the rest. Run it with `step7b_sorosim_finger()`; as a function it doesn't touch your workspace. |

Each script saves a matching `.png` plot in this folder.

## Setup notes
- `sorosim_shims/`: SoRoSim's prebuilt MEX files were compiled with a different MATLAB Coder version and crash. These small wrapper files make SoRoSim call its own MATLAB source instead. Keep this folder; Step 7 adds it to the path automatically.
- `sorosim_shims/CableActuation.m`: replaces SoRoSim's cable-setup window, which saved empty values on this machine. It always creates one straight cable at y = 0, z = −8 mm; edit the two numbers at the top to move it. Delete or rename this file to get the original window back.
- `sorosim_runs/`: SoRoSim writes its solver output here.
- `logs/`: MATLAB-Claude connection logs (safe to ignore).

## Known limits of these models
1. No chamber-wall ballooning (Step 2), so large-bend predictions are too high.
2. Finger springiness (4 Hz, damping 0.15) and valve discharge coefficient are guesses.
3. Valves are modelled as smoothly opening; real solenoids switch on/off with PWM.
4. Grip force uses a lumped base-moment estimate; the real blocked force will differ.
5. The PID gains were found by a coarse grid search; several sit at the edge of the grid.

## What to measure in December (turns this into a calibrated model)
1. Bend angle vs pressure (phone camera + protractor, or an IMU on the tip) → fit the shear modulus in Step 2.
2. Video of a fast inflate/vent → fit the 4 Hz / 0.15 damping guess in Step 4.
3. Tip force vs pressure with a load cell at 2–3 contact angles → check Step 6.
4. Run the Step 5 two-loop controller on an ESP32 with a pressure sensor and a flex sensor or IMU.

Comparing "simple model vs real finger, then calibration" is a small paper's worth of results.

## Next modelling steps
1. Add a point force at the fingertip in SoRoSim ("point wrenches" = 1 when building the linkage) to check the Step 6 grip-force estimate.
2. Switch gravity on to see how much the 36 g finger droops.
3. Once you have December data, fit E in SoRoSim and the shear modulus in Step 2, and compare which model predicts better.
