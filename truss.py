"""
2D pin-jointed truss solver using the direct stiffness method -- the same
underlying technique real structural analysis / FEA software uses (just
simplified to axial-only, pin-jointed members).

Given a truss's geometry (node coordinates), connectivity (which nodes
each member joins, plus its material/section properties), supports
(which node DOFs are fixed), and applied loads, this computes:

  - nodal displacements
  - support reactions
  - the axial force in every member (positive = tension, negative =
    compression)

Method:
  1. For each member, build its 4x4 stiffness matrix in global coordinates
     from its length and orientation.
  2. Assemble all member matrices into one global stiffness matrix K by
     summing contributions at each node's DOF indices.
  3. Partition K into free vs. fixed DOFs, solve K_ff @ u_f = F_f for the
     unknown displacements (fixed DOFs are assumed to have zero
     prescribed displacement -- no support settlement modelled).
  4. Reactions come from the full equilibrium equation K @ u = F + R at
     the fixed DOFs.
  5. Each member's axial force is (EA/L) times its elongation, computed
     directly from the nodal displacements at its two ends.
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


@dataclass
class Member:
    node_i: int
    node_j: int
    E: float  # Young's modulus
    A: float  # cross-sectional area


Node = Tuple[float, float]


def _member_geometry(nodes: List[Node], member: Member):
    xi, yi = nodes[member.node_i]
    xj, yj = nodes[member.node_j]
    dx, dy = xj - xi, yj - yi
    length = (dx**2 + dy**2) ** 0.5
    cx, cy = dx / length, dy / length
    return length, cx, cy


def _member_stiffness_global(member: Member, length, cx, cy) -> np.ndarray:
    k = member.E * member.A / length
    c2, s2, cs = cx * cx, cy * cy, cx * cy
    return k * np.array(
        [
            [c2, cs, -c2, -cs],
            [cs, s2, -cs, -s2],
            [-c2, -cs, c2, cs],
            [-cs, -s2, cs, s2],
        ]
    )


def _dof_indices(node_i, node_j):
    # Each node has 2 DOFs (x, y); global DOF index = 2*node + {0 for x, 1 for y}
    return [2 * node_i, 2 * node_i + 1, 2 * node_j, 2 * node_j + 1]


def solve_truss(
    nodes: List[Node],
    members: List[Member],
    supports: Dict[int, Tuple[bool, bool]],
    loads: Dict[int, Tuple[float, float]],
):
    """
    supports: {node_index: (fix_x, fix_y)} -- True means that DOF is fixed
    loads: {node_index: (Fx, Fy)} -- external applied loads (unspecified nodes get 0)

    Returns:
        displacements: np.ndarray shape (2*n_nodes,) -- [u0, v0, u1, v1, ...]
        reactions: dict {node_index: (Rx, Ry)} for supported nodes
        member_forces: dict {member_index: force} (+ tension, - compression)
    """
    n = len(nodes)
    ndof = 2 * n

    K = np.zeros((ndof, ndof))
    geometry = []
    for member in members:
        length, cx, cy = _member_geometry(nodes, member)
        geometry.append((length, cx, cy))
        k_global = _member_stiffness_global(member, length, cx, cy)
        idx = _dof_indices(member.node_i, member.node_j)
        for a in range(4):
            for b in range(4):
                K[idx[a], idx[b]] += k_global[a, b]

    F = np.zeros(ndof)
    for node_idx, (fx, fy) in loads.items():
        F[2 * node_idx] += fx
        F[2 * node_idx + 1] += fy

    fixed_dofs = set()
    for node_idx, (fix_x, fix_y) in supports.items():
        if fix_x:
            fixed_dofs.add(2 * node_idx)
        if fix_y:
            fixed_dofs.add(2 * node_idx + 1)

    all_dofs = list(range(ndof))
    free_dofs = [d for d in all_dofs if d not in fixed_dofs]
    fixed_dofs = sorted(fixed_dofs)

    K_ff = K[np.ix_(free_dofs, free_dofs)]
    F_f = F[free_dofs]

    if len(free_dofs) > 0:
        u_f = np.linalg.solve(K_ff, F_f)
    else:
        u_f = np.array([])

    displacements = np.zeros(ndof)
    for i, dof in enumerate(free_dofs):
        displacements[dof] = u_f[i]
    # fixed dofs remain 0 (no prescribed settlement)

    full_reaction_forces = K @ displacements  # = F + R everywhere
    reactions = {}
    for node_idx, (fix_x, fix_y) in supports.items():
        rx = full_reaction_forces[2 * node_idx] - F[2 * node_idx] if fix_x else 0.0
        ry = full_reaction_forces[2 * node_idx + 1] - F[2 * node_idx + 1] if fix_y else 0.0
        reactions[node_idx] = (rx, ry)

    member_forces = {}
    for m_idx, member in enumerate(members):
        length, cx, cy = geometry[m_idx]
        idx = _dof_indices(member.node_i, member.node_j)
        u_local = displacements[idx]  # [ui, vi, uj, vj]
        # Elongation = (u_j - u_i) projected onto the member's axis
        elongation = (-cx, -cy, cx, cy) @ u_local
        force = member.E * member.A / length * elongation
        member_forces[m_idx] = force

    return displacements, reactions, member_forces
