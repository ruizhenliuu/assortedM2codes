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

-- Chow ring of the partial blowup of PL at stage m.
-- Works in the matroid Chow ring A(M); returns (A, al, Sk) where
-- Sk contains S_{k,M} for k = 1..min(m+1, r-2), the active stages.
-- m = -1: no blowup (P^{r-1}), Sk is empty.
-- m >= r-2: full wonderful compactification, Sk has all r-2 entries.
partialBlowupChowRing = (M, m) -> (
    r := rank M;
    (A, flatToVar) := makeChowRing M;
    al := alphaClass flatToVar;
    bound := max(0, min(m + 1, r - 2));
    Sk := hashTable for k from 1 to bound list
        k => sum(select(keys flatToVar, F -> rank(M,F) == r-k), F -> flatToVar#F);
    (A, al, Sk))

-- Total Chern class of the tangent bundle of the strict transform of PL,
-- as an element of A(M). Uses truncated recursive formula from
-- arXiv:2510.06609 Corollary 5.5. The blow-up stage is determined by
-- the entries in Sk.
strictTransformChernPoly = (r, A, al, Sk) -> (
    T := new MutableHashTable;
    T#0 = 1_A;
    bound := if #Sk == 0 then 0 else max keys Sk;
    for nn from 1 to bound do (
        ans := (1 + al)^(nn + 1);
        for j from 1 to nn - 1 do
            for i from 0 to nn + 1 - j do
                ans = ans + binomial(nn + 1 - j, i) * T#(j-1) *
                    ((1 + Sk#j) * (1 - Sk#j)^i - 1) *
                    (al - sum(j - 1, kk -> Sk#(kk+1)))^(nn + 1 - i - j);
        T#nn = ans);
    totalChern := (1 + al)^r;
    for j from 1 to bound do
        for i from 0 to r - j do
            totalChern = totalChern + binomial(r - j, i) * T#(j-1) *
                ((1 + Sk#j) * (1 - Sk#j)^i - 1) *
                (al - sum(j - 1, kk -> Sk#(kk+1)))^(r - i - j);
    totalChern)

-- Chern numbers of strict transform of PL at blow-up stage m.
strictTransformChernNumbers = (M, m) -> (
    r := rank M;
    if r < 2 then return {};
    (A, al, Sk) := partialBlowupChowRing(M, m);
    totalChern := strictTransformChernPoly(r, A, al, Sk);
    ck := for k from 0 to r - 1 list (
        b := basis(k, A);
        cs := last coefficients(totalChern, Monomials => b);
        sum(numcols b, j -> (cs_(j,0)) * b_(0,j)));
    degOne := al^(r - 1);
    R2 := ambient A;
    dc := leadCoefficient lift(degOne, R2);
    deg := f -> (leadCoefficient lift(f + degOne/2, R2)) / dc - 1/2;
    apply(partitions(r - 1), p -> (toList p, deg product(toList p, i -> ck#i))))

-- Table of Chern numbers of strict transform at each blow-up stage.
-- Rows indexed by m (-1 to r-2), columns by partitions of r-1.
strictTransformChernTable = M -> (
    r := rank M;
    parts := apply(partitions(r - 1), toList);
    header := prepend("m", apply(parts, p -> concatenate between(",", apply(p, toString))));
    rows := for m from -1 to r - 2 list
        prepend(toString m, apply(strictTransformChernNumbers(M, m), (p, v) -> toString v));
    netList(prepend(header, rows), Alignment => Center, HorizontalSpace => 2))

-- Weak compositions of n into k non-negative parts.
weakCompositions = (k, n) -> (
    if k == 1 then return {{n}};
    result := {};
    for i from 0 to n do
        result = join(result, apply(weakCompositions(k-1, n-i), c -> prepend(i, c)));
    result)

-- Universal Chern number formula for rank r.
-- Works in QQ[z_0..z_{r-1}] where z_0 = alpha, z_k = S_k.
-- c(T_M) = prod_{k=1}^{r-1}(1+z_k) * prod_{i=0}^{r-1}(1+z_0-sum_{j=1}^i z_j).
-- Returns (R, ck, formulas) where ck = Chern class polynomials,
-- formulas = list of (partition, degree r-1 polynomial).
chernNumberFormula = r -> (
    R := QQ[z_0..z_(r-1), Degrees => apply(r, i -> {1})];
    totalChern := product(r-1, k -> 1 + R_(k+1)) *
        product(r, i -> 1 + R_0 - sum(min(i, r-1), j -> R_(j+1)));
    ck := for k from 0 to r-1 list (
        b := basis(k, R);
        cs := last coefficients(totalChern, Monomials => b);
        sum(numcols b, j -> (cs_(j,0)) * b_(0,j)));
    (R, ck, apply(partitions(r-1), p ->
        (toList p, product(toList p, i -> ck#i)))))

-- Mixed intersection numbers deg(alpha^{a_0} * S_1^{a_1} * ... * S_{r-1}^{a_{r-1}})
-- for all weak compositions a of r-1 into r parts.
-- Returns a HashTable mapping exponent lists to degrees.
mixedIntersections = M -> (
    r := rank M;
    (A, flatToVar) := makeChowRing M;
    al := alphaClass flatToVar;
    Sk := for k from 1 to r-1 list
        sum(select(keys flatToVar, F -> rank(M,F) == r-k), F -> flatToVar#F);
    vals := prepend(al, Sk);
    degOne := al^(r-1);
    R2 := ambient A;
    dc := leadCoefficient lift(degOne, R2);
    deg := f -> (leadCoefficient lift(f + degOne/2, R2)) / dc - 1/2;
    hashTable for e in weakCompositions(r, r-1) list
        e => deg product(r, i -> (vals#i)^(e#i)))

-- Chern numbers from the combinatorial formula:
-- c_lambda(M) = sum_a C_{lambda,a} * deg(alpha^{a_0} S_1^{a_1} ... S_{r-1}^{a_{r-1}}).
-- Verifies that the universal formula reproduces matroidChernNumbers.
combinatorialChernNumbers = M -> (
    r := rank M;
    if r < 2 then return {};
    (R, ck, formulas) := chernNumberFormula r;
    mv := mixedIntersections M;
    b := basis(r-1, R);
    apply(formulas, (p, poly) -> (
        cs := last coefficients(poly, Monomials => b);
        cn := sum(numcols b, j -> (
            e := flatten exponents b_(0,j);
            lift(cs_(j,0), QQ) * mv#e));
        (p, cn))))

-- Todd genus via the multiplicative sequence Q(x) = x/(1-e^{-x}).
-- Always 1 for matroid wonderful varieties.
toddGenus = M -> (
    r := rank M; n := r - 1;
    if n == 0 then return 1;
    (A, flatToVar) := makeChowRing M;
    al := alphaClass flatToVar;
    Sk := for k from 1 to r-1 list
        sum(select(keys flatToVar, F -> rank(M,F) == r-k), F -> flatToVar#F);
    qc := {1, 1/2, 1/12, 0, -1/720, 0, 1/30240, 0, -1/1209600};
    Q := x -> sum(min(n+1, #qc), k -> qc#k * x^k);
    toddClass := product(r-1, k -> Q(Sk#k)) *
        product(r, i -> Q(al - sum(min(i,r-1), j -> Sk#j)));
    degOne := al^n; R2 := ambient A;
    dc := leadCoefficient lift(degOne, R2);
    deg := f -> (leadCoefficient lift(f + degOne/2, R2)) / dc - 1/2;
    b := basis(n, A);
    cs := last coefficients(toddClass, Monomials => b);
    deg sum(numcols b, j -> (cs_(j,0)) * b_(0,j)))

-- Hirzebruch chi_y genus. For matroid wonderful varieties (rational),
-- h^{p,q} = 0 for p != q, so chi_y = sum (-1)^p * dim(A^p(M)) * y^p.
-- Specializes to: chi_0 = Todd genus = 1, chi_{-1} = Euler char = c_{r-1}.
chiYGenus = M -> (
    r := rank M; n := r - 1;
    (A, flatToVar) := makeChowRing M;
    R := QQ[y];
    sum(n + 1, p -> (-1)^p * (numcols basis(p, A)) * y^p))
