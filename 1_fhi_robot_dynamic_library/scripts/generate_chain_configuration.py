#!/usr/bin/python3
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
import argparse
import os
import sys

from typing import Dict, List
from dataclasses import dataclass
import yaml
from urdf_parser.paser import Urdf
from urdf_parser.objects import Chain, Link, Joint
from urdf_parser.utils import inertia_matrix_from_vec

# Removes !!python/object: from yaml
yaml.emitter.Emitter.process_tag = lambda *args, **kwargs: None


@dataclass
class Module:
    type: str
    mass: float
    CoM: List[float]
    inertia_matrix: List[List[float]]
    motor_inertia: float
    static_friction: float
    dynamic_friction: float
    gear_ratio: float
    axis: List[float]
    origin: List[float]
    position_limits: List[float]
    speed_limit: float
    effort_limit: float

    @staticmethod
    def from_link_and_joint(link: Link, joint: Joint) -> Module:
        return Module(
            type=joint.type,
            mass=link.mass,
            CoM=link.com.tolist(),
            inertia_matrix=inertia_matrix_from_vec(link.inertias).tolist(),
            motor_inertia=link.motor_inertia,
            static_friction=joint.static_friction,
            dynamic_friction=joint.dynamic_friction,
            gear_ratio=joint.gear_ratio,
            axis=joint.axis.tolist(),
            origin=joint.origin.tolist(),
            position_limits=[joint.limits.lower, joint.limits.upper],
            speed_limit=joint.limits.velocity,
            effort_limit=joint.limits.effort,
        )


@dataclass
class ChainConfiguration:
    modules: Dict[str, Module]

    def from_urdf_string(urdf_string: str) -> ChainConfiguration:
        urdf = Urdf.init_from_string(urdf_string)
        chain = urdf.root_chains[0]
        urdf.merge_fixed_joints(True, "ee_A")
        urdf.print_tree()
        if len(urdf.root_chains) > 1:
            print(
                "Multiple chains found in urdf. "
                "First one has been arbitrarily selected. Check produced yaml file."
            )
        return ChainConfiguration.from_chain(chain)

    def from_chain(chain: Chain) -> ChainConfiguration:
        return ChainConfiguration(
            modules={
                joint.name: Module.from_link_and_joint(link, joint)
                for joint, link in zip(chain.joints, chain.links)
            }
        )

    def _to_yaml(self):
        return yaml.dump(self, sort_keys=False, default_flow_style=None)

    def dump(self, file_name: str):
        with open(file_name, "w") as f:
            f.write(self._to_yaml())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a yaml configuration of a chain from urdf"
    )
    parser.add_argument(
        "--output_file_path",
        default=f"{os.getcwd()}/chain_config.yaml",
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

    with open(args.urdf_file_path, "r") as f:
        chain_configuration = ChainConfiguration.from_urdf_string(f.read())
    chain_configuration.dump(args.output_file_path)
