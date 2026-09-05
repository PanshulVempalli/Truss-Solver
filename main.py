"""
Analyses a 9-node, 15-member Warren truss (a standard footbridge/roof
truss layout) under deck loading, prints the member force table, and
saves a diagram coloured by tension/compression.

Usage:
    python main.py --out warren_truss.png
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from truss import Member, solve_truss

# --- Geometry: a 9-node Warren truss, 8m span, 2m deep ---
NODES = [
    (0.0, 0.0),  # 0 bottom-left  (pin support)
    (2.0, 0.0),  # 1 bottom
    (4.0, 0.0),  # 2 bottom (mid)
    (6.0, 0.0),  # 3 bottom
    (8.0, 0.0),  # 4 bottom-right (roller support)
    (1.0, 2.0),  # 5 top
    (3.0, 2.0),  # 6 top
    (5.0, 2.0),  # 7 top
    (7.0, 2.0),  # 8 top
]

STEEL_E = 200e9  # Pa, structural steel
SECTION_A = 0.002  # m^2, arbitrary uniform section for this example

MEMBER_DEFS = [
    # bottom chord
    (0, 1), (1, 2), (2, 3), (3, 4),
    # top chord
    (5, 6), (6, 7), (7, 8),
    # diagonals (Warren zigzag)
    (0, 5), (5, 1), (1, 6), (6, 2), (2, 7), (7, 3), (3, 8), (8, 4),
]

SUPPORTS = {0: (True, True), 4: (False, True)}  # pin at node 0, roller at node 4

# Deck loads applied at the bottom-chord interior nodes (e.g. footbridge traffic load)
LOADS = {1: (0.0, -12_000.0), 2: (0.0, -12_000.0), 3: (0.0, -12_000.0)}  # N


def build_members():
    return [Member(node_i=i, node_j=j, E=STEEL_E, A=SECTION_A) for i, j in MEMBER_DEFS]


def print_report(reactions, member_forces):
    print("Support reactions:")
    for node, (rx, ry) in reactions.items():
        print(f"  Node {node}: Rx = {rx/1000:8.2f} kN, Ry = {ry/1000:8.2f} kN")

    print("\nMember forces (+ tension, - compression):")
    for idx, (i, j) in enumerate(MEMBER_DEFS):
        force = member_forces[idx]
        kind = "tension" if force > 1e-6 else "compression" if force < -1e-6 else "zero-force"
        print(f"  M{idx:2d}  ({i}-{j}):  {force/1000:8.2f} kN   [{kind}]")


def plot_truss(reactions, member_forces, out_path):
    fig, ax = plt.subplots(figsize=(10, 5))

    max_force = max(abs(f) for f in member_forces.values()) or 1.0

    for idx, (i, j) in enumerate(MEMBER_DEFS):
        xi, yi = NODES[i]
        xj, yj = NODES[j]
        force = member_forces[idx]
        color = "tab:red" if force < 0 else ("tab:blue" if force > 0 else "grey")
        lw = 1.5 + 4.0 * abs(force) / max_force
        ax.plot([xi, xj], [yi, yj], color=color, linewidth=lw, solid_capstyle="round", zorder=1)

    for i, (x, y) in enumerate(NODES):
        ax.plot(x, y, "o", color="black", zorder=2, markersize=5)

    # Mark supports
    ax.plot(*NODES[0], marker="^", color="green", markersize=16, zorder=3, label="Pin support")
    ax.plot(*NODES[4], marker="o", color="orange", markersize=14, zorder=3, label="Roller support")

    ys = [y for _, y in NODES]
    ax.set_ylim(min(ys) - 1.0, max(ys) + 0.5)  # leave room below for load arrows

    # Mark loads
    for node, (fx, fy) in LOADS.items():
        x, y = NODES[node]
        ax.annotate(
            "", xy=(x, y - 0.7), xytext=(x, y - 0.05),
            arrowprops=dict(arrowstyle="->", color="purple", lw=2),
            annotation_clip=False,
        )

    ax.plot([], [], color="tab:blue", lw=3, label="Tension")
    ax.plot([], [], color="tab:red", lw=3, label="Compression")
    ax.plot([], [], color="purple", marker=r"$\downarrow$", linestyle="None", markersize=10, label="Applied load")

    ax.set_title("Warren Truss — Footbridge Analysis (Direct Stiffness Method)")
    ax.set_aspect("equal")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4, fontsize=9)
    ax.set_xlabel("m")
    ax.set_ylabel("m")
    ax.grid(alpha=0.2)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nSaved diagram to {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Warren truss structural analysis")
    parser.add_argument("--out", default="warren_truss.png")
    args = parser.parse_args()

    members = build_members()
    _, reactions, member_forces = solve_truss(NODES, members, SUPPORTS, LOADS)

    print_report(reactions, member_forces)
    plot_truss(reactions, member_forces, args.out)


if __name__ == "__main__":
    main()
