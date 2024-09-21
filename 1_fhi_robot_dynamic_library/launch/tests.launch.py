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
from launch import LaunchDescription
from launch.substitutions import LaunchConfiguration
from launch.actions import DeclareLaunchArgument, OpaqueFunction, RegisterEventHandler
from launch_ros.actions import Node
from fhi_ros2.launch import MoveItConfigsBuilder

from launch.event_handlers import OnProcessExit

TEST_NODES = [
    "dyn_functions_test",
]


def launch_setup(context, *args, **kwargs):
    move_group_node_name = LaunchConfiguration("move_group_node_name").perform(context)
    moveit_remote_config = MoveItConfigsBuilder.get_remote_config(
        context, move_group_name=move_group_node_name
    )

    nodes = [
        Node(
            package="fhi_robot_dynamic_library",
            executable=node_name,
            output="screen",
            prefix=["xterm -e gdb -ex run --args"],
            parameters=[
                moveit_remote_config,
            ],
        )
        for node_name in TEST_NODES
    ]

    events = [
        RegisterEventHandler(
            OnProcessExit(
                target_action=nodes[i],
                on_exit=[nodes[i + 1]],
            )
        )

        for i in range(len(nodes) - 1)
    ]
    return [nodes[0]] + events


def generate_launch_description():
    return LaunchDescription(
        [
            DeclareLaunchArgument(
                "move_group_node_name",
                default_value="move_group",
                description="Name of the 'move_group' node from which info should be published",
            ),
            OpaqueFunction(function=launch_setup),
        ]
    )
