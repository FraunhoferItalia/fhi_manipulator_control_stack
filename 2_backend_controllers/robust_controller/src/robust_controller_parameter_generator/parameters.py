#
# Copyright 2022-2024 Fraunhofer Italia Research
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from __future__ import annotations
import copy
from dataclasses import dataclass
import json
from typing import List
from .objects import Joint, Link
import math
import numpy as np
import scipy
from .utils import maxEigComputation, InertiaMatFromVect


@dataclass
class DynamicParameters:
    m: List[float] = None
    r: List[List[float]] = None
    I: List[List[float]] = None
    Im: List[List[float]] = None
    Fc: List[float] = None
    Fv: List[float] = None
    kr: List[float] = None

    def from_links(links: List[Link]) -> DynamicParameters:
        dyn_par = DynamicParameters()
        dyn_par.m = [link.mass for link in links]
        dyn_par.r = [link.com[:3] for link in links]
        dyn_par.I = [link.inertias for link in links]
        dyn_par.Im = [link.motor_inertia for link in links]
        dyn_par.Fc = [link.static_friction for link in links]
        dyn_par.Fv = [link.dynamic_friction for link in links]
        dyn_par.kr = [link.gear_ratio for link in links]

        return dyn_par

    def __as_dict__(self):
        dict = self.__dict__.copy()
        temp_inertia = np.zeros((3, 3, len(self.m)))
        for i in range(len(self.m)):
            temp_inertia[:, :, i] = InertiaMatFromVect(np.array(self.I[i]))
        dict["I"] = temp_inertia.tolist()

        temp_r = np.zeros((3, len(self.m)))
        for i in range(len(self.m)):
            temp_r[:, i] = np.array(self.r[i])

        dict["r"] = temp_r.tolist()
        return dict


@dataclass
class KinematicParameters:
    xyz: List[List[float]] = None
    rpy: List[List[float]] = None

    def from_joints(joints: List[Joint]) -> KinematicParameters:
        kyn_par = KinematicParameters()
        kyn_par.xyz = [joint.origin[:3] for joint in joints]
        kyn_par.rpy = [joint.origin[3:] for joint in joints]
        return kyn_par

    def __as_dict__(self):
        dict = self.__dict__.copy()
        dict["xyz"] = list(map(list, zip(*self.xyz)))  # tranpose of list
        dict["rpy"] = list(map(list, zip(*self.rpy)))
        return dict


# Non zero initial values may be changed by user at runtime (carefully cit. CN, controller father)
@dataclass
class ControllerParameters:
    robot_parameters: RobotParameters
    Kd: float = 0.0
    Kp: float = 0.0
    Ki: float = 0.0
    PSDmax: float = 0.0
    PSDmin: float = 0.0
    freqLim: float = 25.0 #35.0  # [Hz] frequency limit in Hz
    gainMax: float = 0.0
    kpUnique: float = 1000.0 #40.0
    lambdaK: float = 0.008
    sampleT: float = 0.001
    winSize: int = 128  # window size for FFT
    winPerc: float = 7.0 / 8.0  # window percentage where the counter restarts
    dqThreshold: float = 0.1  # [rad/s], threshold on joint speed for update gain
    B: List[List[float]] = None
    P: List[List[float]] = None
    integral: bool = False

    def __post_init__(self):
        n_dof = self.robot_parameters.nDof
        # basic controller parameter -> PID critically damped
        noise_std_dot = 0.03  # standard deviation of the noise on dq
        sigma_torque = 2 * 0.388579604  # [Nm] maximum acceptable noise on the generated torque
        # compute max Eig of Inertia and avgInertia
        maxEigM, avgDiagM = maxEigComputation(
            self.robot_parameters, 1000
        )  # on 500 monte carlo runs
        lambdaM = 1.4 * maxEigM  # added 40% to be conservative

        if self.integral == False:
            self.Kd = sigma_torque / (noise_std_dot * lambdaM) # lambdaM -> max Eig Inertia
            zeta = 1 # 1 = critically damped
            self.Kp = (self.Kd/(2*zeta))**2

            t_rise = 10 * self.sampleT  # rise time
            omega_max = 3.8897 / t_rise # 3.8897 is the solution of exp(-x)*(1+x) - 0.1
            kd_max = omega_max*2*zeta

            A_1 = np.hstack(
                (np.zeros((n_dof, n_dof)), np.eye(n_dof))
            )
            A_2 = np.hstack(
                (
                    -self.Kp * np.eye(n_dof),
                    -self.Kd * np.eye(n_dof),
                )
            )
            A = np.vstack((A_1, A_2))

            B = np.vstack(
                (np.zeros((n_dof, n_dof)), np.eye(n_dof))
            )

            Q = np.diag(np.hstack((avgDiagM, avgDiagM)))

            P = scipy.linalg.solve_continuous_lyapunov(
                np.transpose(A), -Q
            )  # solution of AX + XA^T = Q

            BP = np.transpose(B) @ P

            self.B = B.tolist()
            self.P = P.tolist()
            self.gainMax = (kd_max - self.Kd) / (BP[-1][-1])


        
        if self.integral == True:
            self.Kd = sigma_torque / (noise_std_dot * lambdaM)  # lambdaM -> max Eig Inertia
            self.Kp = ((self.Kd) ** 2) / 3
            self.Ki = ((self.Kd) ** 3) / 27

            t_rise = 10 * self.sampleT  # rise time
            omega_pole_max = (
                1 / (10 * t_rise) * (2 * math.pi)
            )  # "10 * " added to be more conservative
            kd_max = omega_pole_max * 3

            A_1 = np.hstack(
                (np.zeros((n_dof, n_dof)), np.eye(n_dof), np.zeros((n_dof, n_dof)))
            )
            A_2 = np.hstack(
                (np.zeros((n_dof, n_dof)), np.zeros((n_dof, n_dof)), np.eye(n_dof))
            )
            A_3 = np.hstack(
                (
                    -self.Ki * np.eye(n_dof),
                    -self.Kp * np.eye(n_dof),
                    -self.Kd * np.eye(n_dof),
                )
            )

            A = np.vstack((A_1, A_2, A_3))

            B = np.vstack(
                (np.zeros((n_dof, n_dof)), np.zeros((n_dof, n_dof)), np.eye(n_dof))
            )

            Q = np.diag(np.hstack((avgDiagM, avgDiagM, avgDiagM)))

            P = scipy.linalg.solve_continuous_lyapunov(
                np.transpose(A), -Q
            )  # solution of AX + XA^T = Q

            BP = np.transpose(B) @ P

            self.B = B.tolist()
            self.P = P.tolist()
            self.gainMax = (kd_max - self.Kd) / (BP[-1][-1])

        self.gainMax = 1e04

        sigmaMin = 0.7 #0.4 #0.7  # user defined
        sigmaMax = 0.8 #0.5 #0.8  # user defined

        self.PSDmin = 2 * sigmaMin**2 * self.sampleT
        self.PSDmax = 2 * sigmaMax**2 * self.sampleT

        if (self.freqLim * self.sampleT * self.winSize <= 1):
            raise Exception("freqLim too small with respect to sampleT and winSize") 

    def toDict(self):
        del self.robot_parameters
        return json.loads(
            json.dumps(
                self,
                default=lambda o: o.__dict__,
                sort_keys=False,
                indent=2,
            )
        )


JOINT_TYPE_MAP = {"prismatic": 0, "revolute": 1, "fixed": 2}


class RobotParameters:
    nDof: int
    jType: List[int]
    g: List[float]

    DynPar: DynamicParameters
    KinPar: KinematicParameters
    DynPar_inf: DynamicParameters
    DynPar_sup: DynamicParameters

    def __init__(
        self, joints: List[Joint], links: List[Link], uncertainties: DynamicParameters
    ):
        self.jType = [JOINT_TYPE_MAP[j.type] for j in joints if not j.type == "fixed"]
        self.nDof = len(self.jType)
        self.g = [0, 0, -9.81]  # gravity vector

        self.KinPar = KinematicParameters.from_joints(joints)
        self.DynPar = DynamicParameters.from_links(links)
        self._create_dynamic_param_limits(uncertainties)

    def _create_dynamic_param_limits(self, uncertainties: DynamicParameters):
        self.DynPar_inf = copy.deepcopy(self.DynPar)
        self.DynPar_sup = copy.deepcopy(self.DynPar)
        for key in uncertainties.__dict__:
            item_nominal = getattr(self.DynPar, key)
            item_uncertainties = getattr(uncertainties, key)
            inf, sup = [], []
            for nom, unc in zip(item_nominal, item_uncertainties):
                if isinstance(nom, list):
                    inf.append([n - u * abs(n) for n, u in zip(nom, unc)])
                    sup.append([n + u * abs(n) for n, u in zip(nom, unc)])
                else:
                    inf.append(nom - unc * abs(nom))
                    sup.append(nom + unc * abs(nom))
            setattr(self.DynPar_inf, key, inf)
            setattr(self.DynPar_sup, key, sup)

    def toDict(self):
        return json.loads(
            json.dumps(
                self,
                default=lambda o: o.__as_dict__(),
                sort_keys=False,
                indent=2,
            )
        )

    def __as_dict__(self):
        return self.__dict__.copy()
