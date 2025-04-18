from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from elliptic_algebra import EllipticCurve, GF
from elliptic_algebra.utils import clean_complex

# ----------------------------------------------------------
# CHARACTER SIGNATURES ON THE UNIT CIRCLE 
# ----------------------------------------------------------
# This is a demo of what’s possible with elliptic curve character signatures
# over finite fields. It takes a point on a curve and visualizes its full
# character signature as arrows on the complex unit circle.
#
# Explore, tinker, break things, try weird inputs — this is a playground!
# Use it to develop intuition, spot patterns, or just enjoy the math.
#
# To run: Just execute and check the output PNG image. Change the point,
# curve, or field as you like.

gf = GF(5, [1, 1, 2])
EC = EllipticCurve(gf([1]), gf([2]))
point = EllipticCurve.ECpoint(EC, gf([1]), gf([2]))

# for representation clarity in the plot
def group_keys_by_value(d):
    reverse = defaultdict(list)
    for k, v in d.items():
        reverse[v].append(k)

    result = {}
    for value, keys in reverse.items():
        if len(keys) == 1:
            result[keys[0]] = value
        else:
            result[tuple(keys)] = value
    return result


def plot_signature(g, title="Character Signature"):
    # more manageable size for testing
    print(g.signature)
    fig, ax = plt.subplots(figsize=(40, 40))
    ax.set_aspect('equal')
    ax.set_title(title, fontsize=50)

    # Draw unit circle
    circle = plt.Circle((0, 0), 1, fill=False, color='gray', linestyle='--')
    ax.add_artist(circle)

    # Plot each complex value as an arrow
    sigs = group_keys_by_value(
        dict([(label, clean_complex(x)) for label, x in g.signature]))

    arrow_color = 'blue'
    legend_handles = []

    for i, (label, value) in enumerate(sigs.items(), start=1):
        # Draw the arrow
        ax.arrow(0, 0, value.real, value.imag, head_width=0.05,
                 color=arrow_color, length_includes_head=True)

        # Put number at tip
        ax.text(value.real, value.imag, str(i), fontsize=40,
                ha='center', va='center', color=arrow_color)

        # Legend entry: number and label
        handle = Line2D([0], [0], color=arrow_color,
                        lw=2, label=f"{i}: {label}")
        legend_handles.append(handle)

    # Legend in top right
    ax.legend(handles=legend_handles,
              title='Character functions that map to the same value', loc='upper right')

    # Aesthetic fixes
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    plt.grid(True)

    # Save to file
    plt.savefig(f"{title}.png")
    plt.close()


for g in EC.points_poly[7:8]:
    g = EllipticCurve.ECpoint(EC, g[0], g[1])

    plot_signature(
        g, f"TEST Character Signature of {g} over the curve ({EC.a}, {EC.b}) over GF({gf.p}^{gf.n})")
