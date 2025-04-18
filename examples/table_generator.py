import pandas as pd
from elliptic_algebra import EllipticCurve, GF

# ------------------------------------------------------------------
# SINGULAR CURVE COUNTS PER FIELD ELEMENT 
# ------------------------------------------------------------------
# This script explores how many singular elliptic curves (a, b), each
# element (x, y) of GF(p^n) lies on. For each pair, we compute all
# possible singular curves passing through it, and summarize it in
# a Pandas table — perfect for analysis, intuition, and fun patterns.
# A showcase of proofs inspired by these tables can be found elsewhere.

input = [(5, [1, 1, 2])]
output = []

for p, a in input:
    gf = GF(p, a)
    elements = gf.elements

    if not gf.is_primitive:
        continue

    for x in gf.elements:
        no_curves = []
        for y in elements:
            P = EllipticCurve.FreeECPoint(x, y)
            curves = P.singular_curves()
            no_curves.append(len(curves))
        total_solns = no_curves.count(
            1) + 2 * no_curves.count(2) + 3 * no_curves.count(3)
        no_curves.append(total_solns)
        row = [str(x)] + no_curves
        output.append(row)

        columns = [f"GF({p}^{len(a) - 1})"]+[str(x)
                                             for x in elements] + ["Total"]
        df = pd.DataFrame(output, columns=columns)

        totals = df.iloc[:, 1:].sum(numeric_only=True)
        total_row = ["Total"] + totals.tolist()
        df.loc[len(df)] = total_row

        df.to_csv(f"GF({gf.p}^{gf.n})", index=False)

        print(df)
