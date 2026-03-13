needsPackage "Matroids"
needsPackage "NormalToricVarieties"

-- Constructs the Chow ring A(M) and returns (A, flatToVar).
-- flatToVar: HashTable mapping proper flats to ring variables.
makeChowRing = M -> (
    r := rank M;
    I := idealChowRing M;
    R := ring I;
    A := R/I;
    properFlats := select(flats M, F -> #F > 0 and rank(M,F) < r);
    flatToVar := hashTable apply(#properFlats, i -> properFlats#i => A_i);
    (A, flatToVar))

-- The divisor class alpha = sum of x_F for flats F containing 0.
alphaClass = flatToVar ->
    sum(select(keys flatToVar, F -> member(0, F)), F -> flatToVar#F)

-- Computes the Chern numbers of a matroid M using the formula from
-- arXiv:2510.06609 Theorem 3.5:
--   c(T_M) = prod_{i=1}^{r-1}(1 + S_{i,M}) * prod_{i=0}^{r-1}(1 + alpha - sum_{j=1}^{i} S_{j,M})
-- where S_{k,M} = sum of x_F over rank r-k flats F.
-- Returns a list of pairs (partition of r-1, Chern number).
matroidChernNumbers = M -> (
    r := rank M;
    if r < 2 then return {};
    (A, flatToVar) := makeChowRing M;
    al := alphaClass flatToVar;
    Sk := hashTable for k from 1 to r-1 list
        k => sum(select(keys flatToVar, F -> rank(M,F) == r-k), F -> flatToVar#F);
    totalChern := product(1..r-1, i -> 1 + Sk#i) *
        product(r, i -> 1 + al - sum(min(i, r-1), j -> Sk#(j+1)));
    ck := for k from 0 to r-1 list (
        b := basis(k, A);
        cs := last coefficients(totalChern, Monomials => b);
        sum(numcols b, j -> (cs_(j,0)) * b_(0,j)));
    degOne := al^(r-1);
    R2 := ambient A;
    dc := leadCoefficient lift(degOne, R2);
    deg := f -> (leadCoefficient lift(f + degOne/2, R2)) / dc - 1/2;
    apply(partitions(r-1), p -> (toList p, deg product(toList p, i -> ck#i))))

-- Computes the Chern numbers of a smooth complete toric variety X.
-- For a smooth toric variety, c(TX) = prod_rho (1 + D_rho).
-- Returns a list of pairs (partition, Chern number).
toricChernNumbers = X -> (
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
        prepend(toString m, apply(toricChernNumbers blowupPn(n, m), (p, v) -> toString v));
    netList(prepend(header, rows), Alignment => Center, HorizontalSpace => 2))
