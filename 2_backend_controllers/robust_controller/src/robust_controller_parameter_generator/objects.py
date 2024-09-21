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
from dataclasses import dataclass
from typing import List
from xml.etree.ElementTree import Element
from .utils import (
    recurse_subtree,
    get_subtree,
    rotMat,
    homogMat,
    InertiaMatFromVect,
    Xcross,
    InertiaVectFromMat,
)
import numpy as np


def parse_origin(tree: Element) -> List[float]:
    try:
        result = tree.attrib["xyz"].split(" ")
    except KeyError:
        result = 3 * ["0"]
    
    try:
        result += tree.attrib["rpy"].split(" ")
    except KeyError:
        result += 3 * ["0"]
    return [float(r) for r in result]


@dataclass
class Link:
    xml: Element
    name: str = None
    mass: float = None
    com: List[float] = None
    inertias: List[float] = None
    motor_inertia: float = 0.001
    static_friction: float = 0.0
    dynamic_friction: float = 0.0
    gear_ratio: float = 1.0

    invalid: bool = False

    def __post_init__(self):
        self.name = self.xml.attrib["name"]
        try:
            self.mass = float(
                recurse_subtree(self.xml, ["inertial", "mass"]).attrib["value"]
            )
            self.inertias = self.parse_inertia_tree(
                recurse_subtree(self.xml, ["inertial", "inertia"])
            )
        except AttributeError:
            self.invalid = True
            print(f"skipping link {self.name}")

        try:
            self.com = parse_origin(recurse_subtree(self.xml, ["inertial", "origin"]))
        except AttributeError:
            self.com = [0 for _ in range(6)]

    def parse_inertia_tree(self, tree: Element) -> List[float]:
        return [
            float(tree.attrib["ixx"]),
            float(tree.attrib["iyy"]),
            float(tree.attrib["izz"]),
            float(tree.attrib["ixy"]),
            float(tree.attrib["ixz"]),
            float(tree.attrib["iyz"]),
        ]

    def mergeDyn(self, link: Link, joint: Joint):
        # IMPORTANT: rpy from center of mass position not considered (like if it is always rpz = 000)
        r_com_1 = np.array(self.com[0:3])
        r_com_2 = np.array(link.com[0:3])
        R = rotMat(np.array(joint.origin[3:]))
        H = homogMat(R, np.array(joint.origin[:3]))
        r_com_2_1 = np.transpose(H @ np.transpose(np.hstack((r_com_2, 1))))[0:3]
        r_com_tot = ((r_com_1 * self.mass) + (r_com_2_1 * link.mass)) / (
            self.mass + link.mass
        )
        self.com = [*r_com_tot, 0, 0, 0]  # merged center of mass

        inertia1 = InertiaMatFromVect(np.array(self.inertias))
        inertia2 = InertiaMatFromVect(np.array(link.inertias))

        delta_r_1 = r_com_tot - r_com_1  # difference between r_com_1 and r_com_tot
        delta_r_2 = r_com_tot - r_com_2_1
        inertia1_com = inertia1 + self.mass * np.transpose(Xcross(delta_r_1)) @ Xcross(delta_r_1)
        inertia2_1 = np.transpose(R) @ inertia2 @ R
        inertia2_com = inertia2_1 + link.mass * np.transpose(Xcross(delta_r_2)) @ Xcross(delta_r_2)
        inertia_com_tot = inertia1_com + inertia2_com

        self.inertias = InertiaVectFromMat(inertia_com_tot)  # merged inertias

        self.mass = self.mass + link.mass  # merged mass

        print(f"Link {self.name} merging link {link.name}")


@dataclass
class Joint:
    xml: Element
    name: str = None
    type: str = None
    parent: Link = None
    parent_name: str = None
    child: Link = None
    child_name: str = None
    origin: List[float] = None

    def __post_init__(self):
        self.name = self.xml.attrib["name"]
        self.type = self.xml.attrib["type"]
        self.parent_name = get_subtree(self.xml, "parent").attrib["link"]
        self.child_name = get_subtree(self.xml, "child").attrib["link"]
        try:
            self.origin = parse_origin(get_subtree(self.xml, "origin"))
        except AttributeError:
            self.origin = [0.0] * 6

    def solve(self, links: List[Link]):
        for link in links:
            if link.name == self.parent_name:
                self.parent = link
            elif link.name == self.child_name:
                self.child = link

    def mergeKin(self, joint: Joint):
        R_1 = rotMat(np.array(joint.origin[3:]))
        R_2 = rotMat(np.array(self.origin[3:]))

        H_1 = homogMat(R_1, np.array(joint.origin[:3]))
        H_2 = homogMat(R_2, np.array(self.origin[:3]))

        H_tot = H_1 @ H_2  # compute full rototranslation

        xyz = H_tot[:3, 3]

        roll = np.arctan2(H_tot[2, 1], H_tot[2, 2])
        pitch = np.arctan2(-H_tot[2, 0], np.sqrt(H_tot[2, 1] ** 2 + H_tot[2, 2] ** 2))
        yaw = np.arctan2(H_tot[1, 0], H_tot[0, 0])
        rpy = [roll, pitch, yaw]

        self.origin = [*xyz, *rpy]
