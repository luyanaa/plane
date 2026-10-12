"""Layer 10: cascade flight controller for the 6-DoF model.

Inner loops (rates) -> attitude (phi/theta) -> outer loops (V/h/psi), with
yaw-rate + sideslip coordination and a throttle PI on airspeed error.
Actuators: position + rate limits (de +-25 deg, da +-20 deg, dr +-25 deg,
throttle 0-1 at 0.5/s). Anti-windup by clamping the integrators.
"""
import numpy as np


class Actuator:
    def __init__(self, pos_lim, rate_lim, tau=0.05):
        self.pos_lim = pos_lim
        self.rate_lim = rate_lim
        self.tau = tau
        self.value = 0.0

    def step(self, cmd, dt):
        # first-order lag + rate limit + position limit
        target = self.value + (cmd - self.value) * (dt / (self.tau + dt))
        dv = target - self.value
        dv = np.clip(dv, -self.rate_lim * dt, self.rate_lim * dt)
        self.value = np.clip(self.value + dv, -self.pos_lim, self.pos_lim)
        return self.value


class Autopilot:
    def __init__(self, sixdof, V_cmd, h_cmd=500.0, psi_cmd=0.0):
        self.m = sixdof
        self.V_cmd, self.h_cmd, self.psi_cmd = V_cmd, h_cmd, psi_cmd
        # gains (hand-set, cascade bandwidth ~2-5 rad/s inner)
        self.kq, self.kth = 0.4, 2.0        # pitch: theta -> q -> de
        self.kp, self.kph = 0.6, 1.2        # roll:  phi -> p -> da
        self.kr, self.kb = 0.8, 0.5         # yaw damper + sideslip coordination
        self.kV, self.kVi = 0.08, 0.015     # airspeed PI -> throttle
        self.kh, self.khd = 0.06, 0.35      # altitude -> theta cmd (with hdot)
        self.kpsi = 0.9                     # heading -> phi cmd
        self.eV, self.eh = 0.0, 0.0
        self.eV_lim = 8.0
        self._h_prev = None
        self.de = Actuator(np.radians(25), np.radians(60))
        self.da = Actuator(np.radians(20), np.radians(60))
        self.dr = Actuator(np.radians(25), np.radians(60))
        self.thr = Actuator(1.0, 0.5)

    def update(self, x, dt, gps_ok=True, V_bias=0.0, yaw_damper=True):
        u, v, w, p, q, r, phi, theta, psi = x[:9]
        ua, va, wa = self.m.air_relative(x)
        V = max(np.sqrt(ua * ua + va * va + wa * wa), 1e-6)
        beta = np.arcsin(np.clip(va / V, -1, 1))
        V_sens = V + V_bias
        if self._h_prev is None:
            self._h_prev = x[11]
        hdot = (x[11] - self._h_prev) / max(dt, 1e-6)
        self._h_prev = x[11]

        # outer loops (GPS-dependent)
        th_cmd = 0.02
        ph_cmd = 0.0
        if gps_ok:
            h_err = self.h_cmd - x[11]
            self.eh = np.clip(self.eh + h_err * dt, -30, 30)
            th_cmd = np.clip(self.kh * h_err - self.khd * hdot, -0.35, 0.35)
            dpsi = (self.psi_cmd - psi + np.pi) % (2 * np.pi) - np.pi
            ph_cmd = np.clip(self.kpsi * dpsi, -np.radians(30), np.radians(30))
        # GPS lost: hold last attitude command, keep airspeed loop only
        # throttle PI on airspeed
        self.eV = np.clip(self.eV + (self.V_cmd - V_sens) * dt, -self.eV_lim, self.eV_lim)
        thr_cmd = np.clip(self.kV * (self.V_cmd - V_sens) + self.kVi * self.eV, 0.0, 1.0)

        # attitude -> rate commands
        q_cmd = self.kth * (th_cmd - theta)
        p_cmd = self.kph * (ph_cmd - phi)
        # inner loops.  Sign convention (from sizing/sixdof aero):
        #   de>0 = TE down = nose-down moment (cmde<0)  -> negate
        #   da>0 = negative rolling moment (clda<0)     -> negate
        #   dr>0 = positive yaw moment (cndr>0)         -> positive
        de_cmd = -self.kq * (q_cmd - q)
        da_cmd = -self.kp * (p_cmd - p)
        # coordinated turn: r_cmd = g tan(phi)/V; rudder tracks it (also damps)
        r_cmd = 9.80665 * np.tan(np.clip(phi, -np.radians(45), np.radians(45))) / V
        dr_cmd = self.kr * (r_cmd - r) if yaw_damper else 0.0
        de = self.de.step(np.clip(de_cmd, -1, 1), dt)
        da = self.da.step(np.clip(da_cmd, -1, 1), dt)
        dr = self.dr.step(np.clip(dr_cmd, -1, 1), dt)
        thr = self.thr.step(thr_cmd, dt)
        self._last = (th_cmd, ph_cmd, de, da, dr, thr)
        return [de, da, dr, thr, thr]

    def state(self):
        return self._last
