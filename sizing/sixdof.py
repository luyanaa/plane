"""Layer 9: nonlinear 6-DoF rigid-body dynamics (trim, modes, engine-out).

Same model data as the JSBSim XMLs: VSPAERO aero tables + analytic control
derivatives, engine map + BEMT propeller closure (precomputed full-throttle
thrust grid, throttle-linear), component CG/inertias. JSBSim 1.3.1 (pip)
silently rejects custom aircraft XMLs in this environment (see README), so
the standard rigid-body EOM are integrated here with RK4.

Conventions: body frame at CG, x fwd / y right / z down; Euler 3-2-1.
"""
import numpy as np

from .atmosphere import G, RHO0
from .configs import MISSION
from .geometry import build_geometry
from .mission import Aircraft
from .weight import WeightModel

ELEV_MAX, AIL_MAX, RUD_MAX = np.radians(25.0), np.radians(20.0), np.radians(25.0)
TAU_ELEV, TAU_RUD = 0.5, 0.6
CL_A_H, CL_A_V = 4.5, 3.0


class SixDoF:
    def __init__(self, key, frame, pack):
        self.key, self.frame, self.pack = key, frame, pack
        self.wind = np.zeros(3)          # earth-frame wind, m/s
        self.geo = build_geometry(frame.S, frame.AR, frame.taper, frame.t_over_c,
                                  frame.l_fus, frame.d_fus, n_nacelles=frame.n_nacelles)
        ac = Aircraft(frame, pack, MISSION)
        self.met = ac.metrics(RHO0)
        self.mass = frame.MTOW_kg
        self.W = self.mass * G
        self.cg = self._cg_x(frame, pack, self.geo)
        self.ixx, self.iyy, self.izz = self._inertias(frame, self.geo)
        self._load_aero(key)
        self._load_thrust_grid(ac)

    def set_wind(self, v_earth):
        self.wind = np.array(v_earth, float)

    def air_relative(self, x):
        """Body-frame air-relative velocities from ground body velocities."""
        u, v, w = x[0], x[1], x[2]
        phi, theta, psi = x[6], x[7], x[8]
        cb, sb = np.cos(psi), np.sin(psi)
        st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
        r_eb = np.array([[cb * ct, cb * st * sp - sb * cp, cb * st * cp + sb * sp],
                         [sb * ct, sb * st * sp + cb * cp, sb * st * cp - cb * sp],
                         [-st, ct * sp, ct * cp]])
        w_body = r_eb.T @ self.wind
        return np.array([u - w_body[0], v - w_body[1], w - w_body[2]])

    # ------------------------------------------------------------ setup
    def _cg_x(self, frame, pack, geo):
        _, items = WeightModel().empty_mass(geo, pack.n_eng, pack.P_kw_per_eng)
        x_wing = 0.33 * frame.l_fus
        b_h = float(np.sqrt(geo.S_h * 4.0))
        x_htp = frame.l_fus - 0.35
        x_eng = 0.30 if pack.n_eng == 1 else x_wing - 0.05
        pos = {'wing': x_wing + 0.30 * geo.MAC, 'fuselage': 0.45 * frame.l_fus,
               'HTP': x_htp + 0.30 * geo.S_h / b_h, 'VTP': x_htp + 0.30 * geo.S_h / b_h,
               'landing gear': 0.55 * frame.l_fus, 'engine+mount': x_eng,
               'fuel system': 0.40 * frame.l_fus, 'propellers': x_eng,
               'avionics/FCS': 0.50 * frame.l_fus, 'servos/wiring': 0.55 * frame.l_fus,
               'payload bay': 0.42 * frame.l_fus, 'electrical': 0.55 * frame.l_fus}
        m_e = sum(items.values())
        m_f = frame.MTOW_kg - m_e - MISSION.payload_kg
        return (sum(items[k] * pos[k] for k in items) + (m_f + MISSION.payload_kg)
                * 0.42 * frame.l_fus) / frame.MTOW_kg

    def _inertias(self, frame, geo):
        _, items = WeightModel().empty_mass(geo, self.pack.n_eng, self.pack.P_kw_per_eng)
        ixx = items['wing'] * (geo.b / 2) ** 2 / 3.0 + 2.0
        iyy = 0.10 * frame.MTOW_kg * frame.l_fus ** 2
        return ixx, iyy, ixx + iyy

    def _load_aero(self, key):
        import pandas as pd
        db = pd.read_csv(self._res() / 'vspaero_db.csv')
        d = db[db.prototype == key].sort_values('AoA')
        a = np.radians(d.AoA.values)
        pcl = np.polyfit(a, d.CLtot.values, 1)
        pcm = np.polyfit(a, d.CMytot.values, 1)
        self.alpha = np.concatenate([[-0.35], a, [0.35]])
        self.CL = np.concatenate([[pcl[0] * -0.35 + pcl[1]], d.CLtot.values,
                                  [pcl[0] * 0.35 + pcl[1]]])
        self.Cm = np.concatenate([[pcm[0] * -0.35 + pcm[1]], d.CMytot.values,
                                  [pcm[0] * 0.35 + pcm[1]]])
        cd0 = 0.009 * self.geo.S_wet / self.frame.S
        cd = d.CDtot.values
        self.CD = np.concatenate([[cd[0]], cd, [cd[-1]]]) + (cd0 - cd[0])
        dr = pd.read_csv(self._res() / 'vspaero_derivs.csv')
        r = dr[dr.prototype == key].iloc[0]
        self.cyb, self.clb, self.cnb = r.CY_beta, r.Cl_beta, r.Cn_beta
        sh_s, sv_s = self.geo.S_h / self.frame.S, self.geo.S_v / self.frame.S
        self.cmde = -sh_s * (1.15 / self.geo.MAC) * CL_A_H * TAU_ELEV
        self.clde = sh_s * CL_A_H * TAU_ELEV
        self.cndr = sv_s * (1.15 / self.geo.b) * CL_A_V * TAU_RUD
        self.cydr = sv_s * CL_A_V * TAU_RUD
        self.clda, self.cnda = -0.15, 0.01
        # rate damping (standard tail-volume estimates)
        self.cmq = -2.0 * CL_A_H * sh_s * (1.15 / self.geo.MAC) ** 2
        self.cnr = -2.0 * CL_A_V * sv_s * (1.15 / self.geo.b) ** 2
        self.clp = -0.40

    def _res(self):
        from pathlib import Path
        return Path(__file__).resolve().parents[1] / 'results'

    def _load_thrust_grid(self, ac, n_v=25):
        # momentum-disk linear thrust model (consistent with the takeoff
        # analysis); raw BEMT over-predicts the static region
        p = self.pack
        vmax = max(self.met['V_cruise'] * 1.6, self.met['V_s_TO'] * 2.2)
        self.vgrid = np.linspace(0.5, vmax, n_v)
        t0 = p.static_thrust(RHO0) / p.n_eng
        t_cr = p.eta_p_cruise * p.P_kw_per_eng * 1000.0 / self.met['V_cruise']
        v_cr = self.met['V_cruise']
        t = t0 + (t_cr - t0) * (self.vgrid / v_cr)
        t = np.where(self.vgrid < v_cr, t,
                     p.eta_p_cruise * p.P_kw_per_eng * 1000.0 / self.vgrid)
        self.Tfull = t

    def thrust(self, V, throttle):
        t = np.interp(V, self.vgrid, self.Tfull)
        return max(t * throttle, 0.0)

    # ------------------------------------------------------------ model
    def aero_forces(self, x, ctrl):
        u, v, w, p, q, r, phi, theta, psi = x[:9]
        ua, va, wa = self.air_relative(x)
        V = max(np.sqrt(ua * ua + va * va + wa * wa), 1e-6)
        alpha = np.arctan2(wa, ua)
        beta = np.arcsin(va / V)
        de, da, dr, _, _ = ctrl
        qbar = 0.5 * RHO0 * V * V
        S, b, cbar = self.frame.S, self.geo.b, self.geo.MAC
        cl = np.interp(alpha, self.alpha, self.CL) + self.clde * de
        cd = np.interp(alpha, self.alpha, self.CD)
        cm = np.interp(alpha, self.alpha, self.Cm) + self.cmde * de
        cy = self.cyb * beta + self.cydr * dr
        cll = self.clb * beta + self.clda * da + self.clp * (p * b / (2.0 * V))
        cn = self.cnb * beta + self.cndr * dr + self.cnda * da + self.cnr * (r * b / (2.0 * V))
        cm += self.cmq * (q * cbar / (2.0 * V))
        # wind -> body (standard rotation), small-angle drag/lift decomposition
        ca, sa = np.cos(alpha), np.sin(alpha)
        l_ = qbar * S * cl
        d_ = qbar * S * cd
        y_ = qbar * S * cy
        fx = -d_ * ca + l_ * sa
        fz = -d_ * sa - l_ * ca
        fy = y_
        # moments (stability-axis approximation, linear region)
        mx = qbar * S * b * cll
        my = qbar * S * cbar * cm
        mz = qbar * S * b * cn
        # thrust: nose tractor(s) along +x; twins at y=+-0.75 with yaw moment
        n_eng = self.pack.n_eng
        if n_eng == 1:
            fx += self.thrust(V, ctrl[3])
        else:
            t_l = self.thrust(V, ctrl[3])
            t_r = self.thrust(V, ctrl[4])
            fx += t_l + t_r
            mz += 0.75 * (t_r - t_l)
        return np.array([fx, fy, fz, mx, my, mz])

    def derivs(self, x, ctrl):
        u, v, w, p, q, r, phi, theta, psi, xe, ye, h = x
        F = self.aero_forces(x, ctrl)
        # gravity in body axes (z down)
        st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
        F[0] += -self.W * st
        F[1] += self.W * sp * ct
        F[2] += self.W * cp * ct
        ud, vd, wd = F[0] / self.mass + r * v - q * w, \
                     F[1] / self.mass - r * u + p * w, \
                     F[2] / self.mass + q * u - p * v
        pd = (F[3] + (self.iyy - self.izz) * q * r) / self.ixx
        qd = (F[4] + (self.izz - self.ixx) * p * r) / self.iyy
        rd = (F[5] + (self.ixx - self.iyy) * p * q) / self.izz
        phid = p + q * np.tan(theta) * np.sin(phi) + r * np.tan(theta) * np.cos(phi)
        thetad = q * np.cos(phi) - r * np.sin(phi)
        psid = (q * np.sin(phi) + r * np.cos(phi)) / max(np.cos(theta), 1e-6)
        # wind->earth velocity (z down -> h up)
        Vb = np.array([u, v, w])
        cb, sb = np.cos(psi), np.sin(psi)
        R_eb = np.array([[cb * ct, cb * st * sp - sb * cp, cb * st * cp + sb * sp],
                         [sb * ct, sb * st * sp + cb * cp, sb * st * cp - cb * sp],
                         [-st, ct * sp, ct * cp]])
        Ve = R_eb @ Vb
        return np.array([ud, vd, wd, pd, qd, rd, phid, thetad, psid,
                         Ve[0], Ve[1], -Ve[2]])

    def integrate(self, x0, ctrl, t_end, dt=0.02):
        x = np.array(x0, float)
        t = 0.0
        while t < t_end:
            k1 = self.derivs(x, ctrl)
            k2 = self.derivs(x + 0.5 * dt * k1, ctrl)
            k3 = self.derivs(x + 0.5 * dt * k2, ctrl)
            k4 = self.derivs(x + dt * k3, ctrl)
            x += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
        return x

    # ------------------------------------------------------------ trim
    def trim_level(self, V, throttle0=0.3, tol=1e-8):
        """Solve (alpha, de, throttle) for steady level flight at V (theta=alpha)."""
        def residuals(z):
            alpha, de, thr = z
            x = np.array([V * np.cos(alpha), 0.0, V * np.sin(alpha),
                          0, 0, 0, 0, alpha, 0, 0, 0, 500.0])
            d = self.derivs(x, [de, 0.0, 0.0, thr, thr])
            return np.array([d[0], d[2], d[4]])   # u-dot, w-dot, q-dot
        z = np.array([0.05, 0.0, throttle0])
        for _ in range(80):
            r = residuals(z)
            if np.max(np.abs(r)) < tol:
                break
            j = np.zeros((3, 3))
            for i in range(3):
                zp = z.copy()
                zp[i] += 1e-6
                j[:, i] = (residuals(zp) - r) / 1e-6
            dz = np.linalg.solve(j, -r)
            z += np.clip(dz, -0.3, 0.3)
        return np.array([z[0], z[1], z[0], z[2]]), residuals(z)

    def trim_controls(self, V, throttle0=0.3):
        z, r = self.trim_level(V, throttle0)
        return {'V': V, 'alpha': z[0], 'elevator': z[1], 'theta': z[2],
                'throttle': z[3], 'residual': np.max(np.abs(r))}
