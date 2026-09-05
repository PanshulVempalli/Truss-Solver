import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from truss import Member, solve_truss


def test_simple_triangle_matches_hand_calculation():
    """
    Classic statically-determinate 3-bar truss:

        C (3,4)
       /|\\
      / | \\
     A--+---B    A=(0,0) pin, B=(6,0) roller (vertical only)

    Members: AB (base), AC, BC. Downward load P at apex C.

    By symmetry and hand calculation (method of joints):
      F_AC = F_BC = -0.625*P   (compression)
      F_AB = +0.375*P          (tension)
      Reactions: Ray = Rby = P/2, Rax = 0

    This is independent of E and A for a statically determinate truss --
    the stiffness method should reproduce these exact numbers regardless
    of what stiffness values are used, which this test also checks by
    running it twice with different E*A.
    """
    nodes = [(0.0, 0.0), (6.0, 0.0), (3.0, 4.0)]  # A=0, B=1, C=2
    P = 100.0

    supports = {0: (True, True), 1: (False, True)}  # A: pin, B: roller (vertical)
    loads = {2: (0.0, -P)}

    for E, A in [(200e9, 0.001), (70e9, 0.0005)]:  # steel-ish, then aluminium-ish
        members = [
            Member(node_i=0, node_j=1, E=E, A=A),  # AB
            Member(node_i=0, node_j=2, E=E, A=A),  # AC
            Member(node_i=1, node_j=2, E=E, A=A),  # BC
        ]
        _, reactions, member_forces = solve_truss(nodes, members, supports, loads)

        assert abs(member_forces[0] - 0.375 * P) < 1e-6, "AB should be 0.375P tension"
        assert abs(member_forces[1] - (-0.625 * P)) < 1e-6, "AC should be 0.625P compression"
        assert abs(member_forces[2] - (-0.625 * P)) < 1e-6, "BC should be 0.625P compression"

        rax, ray = reactions[0]
        rbx, rby = reactions[1]
        assert abs(rax) < 1e-6, "pin at A should have zero horizontal reaction by symmetry"
        assert abs(ray - P / 2) < 1e-6
        assert abs(rby - P / 2) < 1e-6


def test_equilibrium_is_satisfied():
    """Sanity check independent of the analytical solution: sum of all
    reactions + applied loads must equal zero (global equilibrium),
    for an arbitrary (still determinate) truss and load case."""
    nodes = [(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]
    members = [
        Member(0, 1, E=200e9, A=0.002),
        Member(1, 2, E=200e9, A=0.002),
        Member(2, 3, E=200e9, A=0.002),
        Member(3, 0, E=200e9, A=0.002),
        Member(0, 2, E=200e9, A=0.002),  # diagonal bracing makes it determinate
    ]
    supports = {0: (True, True), 1: (False, True)}
    loads = {3: (50.0, -80.0), 2: (0.0, -40.0)}

    _, reactions, _ = solve_truss(nodes, members, supports, loads)

    total_fx = sum(fx for fx, _ in loads.values()) + sum(r[0] for r in reactions.values())
    total_fy = sum(fy for _, fy in loads.values()) + sum(r[1] for r in reactions.values())

    assert abs(total_fx) < 1e-6, "horizontal forces must balance"
    assert abs(total_fy) < 1e-6, "vertical forces must balance"


def test_zero_load_gives_zero_forces():
    """With no applied load, every member force and reaction should be zero."""
    nodes = [(0.0, 0.0), (6.0, 0.0), (3.0, 4.0)]
    members = [
        Member(0, 1, E=200e9, A=0.001),
        Member(0, 2, E=200e9, A=0.001),
        Member(1, 2, E=200e9, A=0.001),
    ]
    supports = {0: (True, True), 1: (False, True)}
    loads = {}

    _, reactions, member_forces = solve_truss(nodes, members, supports, loads)

    for f in member_forces.values():
        assert abs(f) < 1e-9
    for rx, ry in reactions.values():
        assert abs(rx) < 1e-9 and abs(ry) < 1e-9


def test_tension_and_compression_sign_convention():
    """A single horizontal member pulled apart at both ends by equal and
    opposite forces should read as pure tension (positive force)."""
    nodes = [(0.0, 0.0), (2.0, 0.0)]
    members = [Member(0, 1, E=200e9, A=0.001)]
    supports = {0: (True, True), 1: (False, True)}
    loads = {1: (10.0, 0.0)}

    _, _, member_forces = solve_truss(nodes, members, supports, loads)
    assert member_forces[0] > 0, "pulling a member apart should register as tension (positive)"
