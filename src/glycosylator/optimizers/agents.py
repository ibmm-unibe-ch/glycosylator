"""
Intelligent agents for optimizing molecules
"""

from typing import Union
import numpy as np
import glycosylator.optimizers.environments as env


class MemoryRotator:
    """
    This class is used to run an agent until convergence
    with a tolerance for suboptimality. It will sample random rotations
    around edges and accept them if they improve the score.

    Parameters
    ----------
    steps : int, optional
        The maximal number of steps to run the agent for, by default 10000

    tolerance : int, optional
        The number of steps to run the agent for to accept convergence, by default 50
        This is also the number of steps the agent will accept a suboptimal score.
    """

    def __init__(self, steps: int = 10000, tolerance: int = 100):
        self.steps = steps
        self.tolerance = tolerance

    def sample(self):
        """
        Sample a random action

        Returns
        -------
        tuple
            The action to take
        """
        rewards = np.ones(len(self._bond_choices))
        # rewards = self._bond_rewards - np.min(self._bond_rewards) + 1
        rewards = rewards / np.sum(rewards)
        bond = np.random.choice(self._bond_choices, p=rewards)

        rewards = self._angle_rewards[2] - np.min(self._angle_rewards[2]) + 1
        rewards = rewards / np.sum(rewards)
        self.adx = np.random.choice(self._angle_choices, p=rewards)

        angle = np.random.uniform(*self._angle_rewards[:2, self.adx])

        return bond, angle

    def update(self, bond, reward):
        """
        Update the mean rewards for the given bond and angle

        Parameters
        ----------
        bond : int
            The bond to update
        angle : float
            The angle to update
        reward : float
            The reward to update the mean with
        """
        self._angle_rewards[2, self.adx] += reward
        self._bond_rewards[bond] += reward

    def run(
        self,
        environment,
        noise: float = 0.01,
        render: Union[int, float] = None,
        render_kws=None,
    ):
        """
        Run the agent until convergence

        Parameters
        ----------
        environment : Environment
            The environment to run the agent in.
        noise : float, optional
            The probability of taking a random action, by default 0.05.

        render : int or float
            The interval at which to render the environment.
            This can be either a fixed number of steps or a fraction of the total number of steps
            By default rendering is disabled.
        render_kws : dict
            The keyword arguments to pass to the render method of the environment
        Returns
        -------
        Environment
            The converged environment
        policy : list
            The policy that was used to converge the environment
        """

        state, reward, *_ = environment.step((0, 0))

        worst_reward = np.inf

        tolerance_counter = 0
        policy = {i: 0 for i in range(len(environment.bonds))}
        temporary_policy = dict(policy)

        last_favourable_state = state
        last_favourable_reward = reward
        reward_history = np.full(20, worst_reward)
        reward_history[0] = reward

        self.setup_sampling_space(environment)

        if isinstance(render, float):
            render = int(render * self.steps)
        if not render_kws:
            render_kws = {}

        # Run the agent until convergence
        for i in range(0, int(self.steps)):

            # Take a step
            action = self.sample()
            new_state, reward, done, info = environment.step(action)

            if reward > last_favourable_reward:

                self.update(action[0], reward * 0.5)

                policy[action[0]] += action[1]

                last_favourable_state = deepcopy(new_state)
                last_favourable_reward = reward

                reward_history = np.roll(reward_history, 1)
                reward_history[0] = reward

                tolerance_counter = 0

                for k in temporary_policy:
                    policy[k] += temporary_policy[k]
                    temporary_policy[k] = 0

            elif reward > worst_reward:

                self.update(action[0], reward * 0.05)
                temporary_policy[action[0]] += action[1]
                tolerance_counter += 1

            else:
                self.update(action[0], -reward * 0.02)
                worst_reward = reward
                tolerance_counter += 1

            if tolerance_counter > self.tolerance:
                tolerance_counter = 0
                environment.set_state(deepcopy(last_favourable_state))
                for k in temporary_policy:
                    temporary_policy[k] = 0

            if (
                np.sum(reward_history / np.max(reward_history) > 0.95)
                > 0.6 * self.tolerance
            ):
                break

        environment.set_state(last_favourable_state)
        return environment, policy, last_favourable_reward

    def setup_sampling_space(self, environment):
        self._bond_choices = np.arange(0, len(environment.bonds))
        self._bond_rewards = np.zeros(len(environment.bonds))

        angle_linspace = np.linspace(0, np.pi, 101)
        self._angle_rewards = np.array(
            [
                [angle_linspace[i] for i in range(len(angle_linspace) - 1)],
                [angle_linspace[i + 1] for i in range(len(angle_linspace) - 1)],
                np.zeros(len(angle_linspace) - 1),
            ]
        )
        self._angle_choices = np.arange(0, 100)


if __name__ == "__main__":

    import glycosylator as gl
    from time import time
    import matplotlib.pyplot as plt

    mol = gl.Molecule.from_pdb(
        "/Users/noahhk/GIT/glycosylator/support/examples/man9.pdb"
    )
    mol.infer_bonds(restrict_residues=False)

    from copy import deepcopy

    mol2 = deepcopy(mol)

    connections = mol.get_residue_connections()
    _cons = mol.get_residue_connections()

    for i in range(len(connections)):
        random_bond = _cons.pop()
        for i in range(5):
            random_angle = np.random.uniform(0.2, np.pi)
            mol.rotate_around_bond(
                *random_bond, np.degrees(random_angle), descendants_only=True
            )

    v = gl.utils.visual.MoleculeViewer3D(mol)
    v.draw_edges(mol2.bonds, color="cyan", linewidth=2)

    graph = mol.make_atom_graph() #mol.make_residue_graph(detailed=True)

    env = env.SingleBondRotatron(graph, connections)
    # env.mask_effector_coords(
    #     np.array([i in graph.residues for i in graph.nodes], dtype=bool)
    # )

    initial_reward = env._compute_reward(env._coords)

    print("----------" * 4)
    conv = MemoryRotator(steps=1000, tolerance=10)
    start = time()
    env, policy, reward = conv.run(env)

    print(
        "Initial reward:",
        initial_reward,
        "Final reward:",
        env._compute_reward(env._coords),
    )
    print(time() - start)

    # apply the policy
    for k, angle in policy.items():
        bond = env.bonds[k]
        angle = np.degrees(angle)
        mol.rotate_around_bond(*bond, angle, descendants_only=True)

    v.draw_edges(mol.bonds, color="magenta")
    v.show()

    v = gl.utils.visual.MoleculeViewer3D(mol)
    v.show()