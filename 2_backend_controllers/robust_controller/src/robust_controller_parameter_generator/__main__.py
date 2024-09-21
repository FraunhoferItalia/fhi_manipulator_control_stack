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
import json
import os
import sys
from typing import List
import xml.etree.ElementTree as ET
from .objects import Link, Joint
from .parameters import (
    RobotParameters,
    DynamicParameters,
    ControllerParameters,
)


def main(input_urdf: str, output_path: str):
    links: List[Link] = []
    joints: List[Joint] = []

    with open(input_urdf, "r") as f:
        urdf_root = ET.fromstring(f.read())

    # Parse urdf creating a first class description
    for child in urdf_root:
        if child.tag == "link":
            link = Link(child)
            # Ignore dynamic-less links
            if not link.invalid:
                links.append(link)
        elif child.tag == "joint":
            joints.append(Joint(child))

    # Solve for joints and links connections based on names
    for joint in reversed(joints):
        joint.solve(links)
        # Invalidate joints without a child or a parent
        if joint.child is None or joint.parent is None:
            joints.pop(joints.index(joint))

    # Remove fixed joints excluding the last one.
    fixed_joints: List[Joint] = [joint for joint in joints if joint.type == "fixed"]
    if fixed_joints:
        for fixed_joint in reversed(fixed_joints):
            fixed_joint.parent.mergeDyn(fixed_joint.child, fixed_joint)
            # mergeKin currently not implemented because with Pino rotor and stator have null xyz and rpz
            links.pop(links.index(fixed_joint.child))
            joints.pop(joints.index(fixed_joint))

        for joint in reversed(joints):
            if joint.type == "revolute":
                prev_fixed_joint = None
                for j in fixed_joints:
                    if j.child == joint.parent:
                        prev_fixed_joint = j
                if not prev_fixed_joint is None:
                    joint.mergeKin(prev_fixed_joint)

        joints.append(fixed_joints[-1])

    print([link.name for link in links])
    print([joint.name for joint in joints])

    N_JOINTS = len(joints)
    uncertainties = DynamicParameters(
        m=[0.1] * N_JOINTS,
        r=[[0.1, 0.1, 0.1]] * N_JOINTS,
        I=[[0.1, 0.1, 0.1, 0.1, 0.1, 0.1]] * N_JOINTS,
        Im=[0.0] * N_JOINTS,
        Fc=[0.0] * N_JOINTS,
        Fv=[0.0] * N_JOINTS,
        kr=[0.0] * N_JOINTS,
    )

    for l in links:
        l.static_friction = 1.3526947430468768
        l.dynamic_friction = 2.0594395737675875

    robot_parameters = RobotParameters(
        joints=joints,
        links=links,
        uncertainties=uncertainties,
    )

    robot_param_dict = robot_parameters.toDict()

    controller_parameters = ControllerParameters(robot_parameters)
    controller_param_dict = controller_parameters.toDict()

    with open(os.path.join(output_path), "w") as f:
        f.write(
            json.dumps(
                {**robot_param_dict, **controller_param_dict},
                sort_keys=False,
                indent=2,
            )
        )

import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create the json file for the robust controller"
    )
    parser.add_argument(
        "--output_file_path",
        default="",
        help="Absolute path of the generated chain configuration (default: <current working directory>/chain_config.yaml)",
    )
    parser.add_argument(
        "--urdf_file_path",
        default="",
        help='Absolute path of robot description urdf (default: "").\nNote: it must be already parsed if xacro',
    )

    args = parser.parse_args()
    if args.urdf_file_path == "":
        print("Cannot operate without a urdf file. Exiting.")
        sys.exit(1)
    if args.output_file_path == "":
        print("--output_file_path not specified.")
        print(
            f"Creating file in default location: {os.getcwd()}/robust_controller_config.yaml"
        )

    main(args.urdf_file_path, args.output_file_path)
