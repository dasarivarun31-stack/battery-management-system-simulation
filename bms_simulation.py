"""Battery Management System (BMS) Simulation
- 1-RC Thevenin cell model + lumped thermal model
- SOC: Coulomb counting vs Extended Kalman Filter (EKF)
- SOH: capacity-based estimate from EKF SOC
- Thermal management: fan cooling + protection cut-offs
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(1)
dt = 1.0
# ---------- OCV curve ----------
soc_pts = np.array([0, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
ocv_pts = np.array([3.00, 3.35, 3.48, 3.58, 3.64, 3.69, 3.75, 3.84, 3.93, 4.02, 4.11, 4.20])
ocv = lambda s: np.interp(s, soc_pts, ocv_pts)
docv = lambda s: (ocv(min(s + .005, 1)) - ocv(max(s - .005, 0))) / (min(s + .005, 1) - max(s - .005, 0))

# ---------- Cell / pack parameters ----------
Q_NOM = 10.0           # Ah rated
SOH_TRUE = 0.92        # aged cell
Q_TRUE = Q_NOM * SOH_TRUE
R0, R1, C1 = 0.012, 0.008, 2500.0
Cth, Rth_nat, Rth_fan = 250.0, 4.0, 1.4   # J/K, K/W
T_amb = 30.0
# ---------- BMS limits ----------
V_MAX, V_MIN, T_FAN_ON, T_FAN_OFF, T_CUT, I_MAX = 4.20, 3.00, 35.0, 32.0, 55.0, 35.0

# ---------- Current profile (+ = discharge) ----------
N = 4200
t = np.arange(N) * dt
I_cmd = np.zeros(N)
I_cmd[:2200] = 11 + 8 * np.sin(2 * np.pi * t[:2200] / 240) + 3 * np.sin(2 * np.pi * t[:2200] / 37)
I_cmd[2200:2500] = 0
I_cmd[2500:] = -10.0   # CC charge

# ---------- State ----------
soc_true, v1, T = 0.90, 0.0, T_amb
soc_cc = 0.60           # coulomb counter (wrong initial guess)
x = np.array([0.60, 0.0])  # EKF state
P = np.diag([0.1, 1e-3]); Qn = np.diag([1e-7, 1e-6]); Rn = 4e-4
fan = False; fault = ""
log = {k: [] for k in "I V Vm soc soc_cc soc_ekf T fan".split()}
ah_dis = 0.0; ekf_s0 = None

for k in range(N):
    I = I_cmd[k]
    # ---- BMS protection ----
    Vlast = log["V"][-1] if log["V"] else ocv(soc_true)
    fault = ""
    if T >= T_CUT: fault = "OVERTEMP"
    if I > 0 and Vlast <= V_MIN: fault = "UNDERVOLT"
    if I < 0 and Vlast >= V_MAX: fault = "OVERVOLT(charge done)"
    if abs(I) > I_MAX: I = np.sign(I) * I_MAX
    if fault: I = 0.0
    # ---- Thermal management (hysteresis) ----
    if T > T_FAN_ON: fan = True
    if T < T_FAN_OFF: fan = False
    Rth = Rth_fan if fan else Rth_nat
    # ---- Plant ----
    soc_true -= I * dt / (3600 * Q_TRUE)
    e = np.exp(-dt / (R1 * C1))
    v1 = v1 * e + R1 * (1 - e) * I
    V = ocv(soc_true) - v1 - R0 * I
    Qgen = I**2 * (R0 + R1)
    T += dt / Cth * (Qgen - (T - T_amb) / Rth)
    # ---- Sensors ----
    Vm = V + rng.normal(0, 0.003); Im = I + 0.15 + rng.normal(0, 0.05)
    # ---- SOC: Coulomb counting ----
    soc_cc -= Im * dt / (3600 * Q_NOM)
    # ---- SOC: EKF ----
    x_p = np.array([x[0] - Im * dt / (3600 * Q_NOM), x[1] * e + R1 * (1 - e) * Im])
    F = np.array([[1, 0], [0, e]])
    P = F @ P @ F.T + Qn
    H = np.array([[docv(x_p[0]), -1.0]])
    y = Vm - (ocv(x_p[0]) - x_p[1] - R0 * Im)
    S = H @ P @ H.T + Rn; K = P @ H.T / S
    x = x_p + (K * y).ravel(); x[0] = np.clip(x[0], 0, 1)
    P = (np.eye(2) - K @ H) @ P
    # ---- SOH bookkeeping (discharge segment, after EKF converged) ----
    if k == 400: ekf_s0 = x[0]; ah_dis = 0.0
    if 400 < k < 2200 and Im > 0: ah_dis += Im * dt / 3600
    if k == 2199: ekf_s1 = x[0]
    for n, v in zip(log, [I, V, Vm, soc_true, soc_cc, x[0], T, fan]): log[n].append(v)

Q_est = ah_dis / (ekf_s0 - ekf_s1)
SOH_est = Q_est / Q_NOM
L = {k: np.array(v) for k, v in log.items()}
rmse = lambda a: np.sqrt(np.mean((a[400:] - L["soc"][400:])**2)) * 100

print("=== BMS SIMULATION RESULTS ===")
print(f"SOC RMSE  Coulomb counting : {rmse(L['soc_cc']):.2f} %")
print(f"SOC RMSE  EKF (t>400 s)    : {rmse(L['soc_ekf']):.2f} %")
print(f"SOH true / estimated       : {SOH_TRUE*100:.1f} % / {SOH_est*100:.1f} %")
print(f"Peak temperature           : {L['T'].max():.1f} C   (fan ON {L['fan'].mean()*100:.0f}% of time)")
print(f"Voltage range              : {L['V'].min():.2f} - {L['V'].max():.2f} V")
print(f"Final BMS state            : {fault or "NORMAL"}")

# ---------- Plots ----------
fig, ax = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
ax[0].plot(t, L["I"], "k"); ax[0].set_ylabel("Current (A)"); ax[0].set_title("Battery Management System Simulation")
ax[1].plot(t, L["Vm"], lw=.5, c="gray", label="Measured"); ax[1].plot(t, L["V"], "b", lw=1, label="True")
ax[1].axhline(V_MIN, c="r", ls="--", lw=.8); ax[1].axhline(V_MAX, c="r", ls="--", lw=.8); ax[1].set_ylabel("Voltage (V)"); ax[1].legend()
ax[2].plot(t, L["soc"]*100, "k", label="True"); ax[2].plot(t, L["soc_cc"]*100, "r--", label="Coulomb counting")
ax[2].plot(t, L["soc_ekf"]*100, "g", label="EKF"); ax[2].set_ylabel("SOC (%)"); ax[2].legend()
ax[3].plot(t, L["T"], "darkorange", label="Cell temp"); ax[3].fill_between(t, 25, 25 + L["fan"]*2, color="c", alpha=.6, label="Fan ON")
ax[3].axhline(T_FAN_ON, c="b", ls=":", lw=.8); ax[3].axhline(T_CUT, c="r", ls="--", lw=.8)
ax[3].set_ylabel("Temp (°C)"); ax[3].set_xlabel("Time (s)"); ax[3].legend()
for a in ax: a.grid(alpha=.3)
plt.tight_layout(); plt.savefig("/mnt/user-data/outputs/bms_results.png", dpi=130)
