needsPackage "NormalToricVarieties"

-- Computes the Chern numbers of a smooth complete toric variety X.
-- For a smooth toric variety, c(TX) = prod_rho (1 + D_rho).
-- Returns a list of pairs (partition, Chern number).
chernNumbers = X -> (
    (n, A, nR) := (dim X, intersectionRing X, #rays X);
    ck := toList fold(
        (c, i) -> apply(#c, k -> if k == 0 then c#0 else c#k + c#(k-1) * A_i),
        apply(n + 1, k -> if k == 0 then 1_A else 0_A),
        0..nR-1);
    pt := product apply(toList(max X)#0, i -> A_i);
    deg := f -> (leadCoefficient lift(f + pt/2, ambient A)) /
                (leadCoefficient lift(pt, ambient A)) - 1/2;
    apply(partitions n, p -> (toList p, deg product(toList p, i -> ck#i))))

-- Iterated blowup of P^n at all coordinate subspaces of dimension 0 through m.
-- When m = -1, this gives projective space P^n.
-- When m = n-1, this gives the permutohedral variety.
blowupPn = (n, m) -> (
    X := toricProjectiveSpace n;
    ground := toList(0..n);
    for d from 0 to min(m, n-2) do
        for S in subsets(ground, d + 1) do
            X = toricBlowup(sort select(ground, i -> not member(i, S)), X);
    X)

permutohedralVariety = n -> blowupPn(n, n-1)

-- Table of Chern numbers of blowupPn(n, m) for m = -1, ..., n-1.
-- Rows are indexed by m, columns by partitions of n.
chernTable = n -> (
    parts := apply(partitions n, toList);
    header := prepend("m", apply(parts, p -> concatenate between(",", apply(p, toString))));
    rows := for m from -1 to n-1 list
        prepend(toString m, apply(chernNumbers blowupPn(n, m), (p, v) -> toString v));
    netList(prepend(header, rows), Alignment => Center, HorizontalSpace => 2))
