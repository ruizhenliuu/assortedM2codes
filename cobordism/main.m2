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

-- Choose a distinguished non-loop to define alpha; any element gives the same class.
distinguishedElement = M -> first sort toList(M.groundSet)

-- The divisor class alpha = sum of x_F for flats F containing the distinguished element.
alphaClass = (M, flatToVar) ->
    sum(select(keys flatToVar, F -> member(distinguishedElement M, F)), F -> flatToVar#F)

-- Computes the Chern numbers of a matroid M using the formula from
-- arXiv:2510.06609 Theorem 3.4:
--   c(T_M) = prod_{i=1}^{r-1}(1 + S_{i,M}) * prod_{i=0}^{r-1}(1 + alpha - sum_{j=1}^{i} S_{j,M})
-- where S_{k,M} = sum of x_F over rank r-k flats F.
-- Returns a list of pairs (partition of r-1, Chern number).
matroidChernNumbers = M -> (
    r := rank M;
    if r < 2 then return {};
    d := r - 1;
    (A, flatToVar) := makeChowRing M;
    al := alphaClass(M, flatToVar);
    Sk := hashTable for k from 1 to d list
        k => sum(select(keys flatToVar, F -> rank(M,F) == r-k), F -> flatToVar#F);
    W := join(
        apply(toList(1..d), i -> Sk#i),
        apply(toList(0..d), j -> if j == 0 then al else al - sum(1..j, i -> Sk#i)));
    E := new MutableList from apply(d + 1, k -> if k == 0 then 1_A else 0_A);
    scan(W, w -> scan(reverse toList(1..d), k -> E#k = E#k + w * E#(k-1)));
    ck := toList E;
    degOne := al^d;
    R2 := ambient A;
    dc := leadCoefficient lift(degOne, R2);
    deg := f -> (leadCoefficient lift(f + degOne/2, R2)) / dc - 1/2;
    apply(partitions d, p -> (toList p, deg product(toList p, i -> ck#i))))

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
    al := alphaClass(M, flatToVar);
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
-- Uses DP to compute elementary symmetric polynomials E_k = e_k[W_r] where
-- W_r = (z_1,...,z_{r-1}, z_0, z_0-z_1, ..., z_0-z_1-...-z_{r-1}).
-- Returns (R, ck, formulas) where ck = Chern class polynomials,
-- formulas = list of (partition, degree r-1 polynomial).
chernNumberFormula = r -> (
    R := QQ[z_0..z_(r-1), Degrees => apply(r, i -> {1})];
    d := r - 1;
    W := join(
        apply(toList(1..d), i -> R_i),
        apply(toList(0..d), j -> if j == 0 then R_0 else R_0 - sum(1..j, i -> R_i)));
    E := new MutableList from apply(d + 1, k -> if k == 0 then 1_R else 0_R);
    scan(W, w -> scan(reverse toList(1..d), k -> E#k = E#k + w * E#(k-1)));
    ck := toList E;
    (R, ck, apply(partitions d, p ->
        (toList p, product(toList p, i -> ck#i)))))

-- Mixed intersection numbers deg(alpha^{a_0} * S_1^{a_1} * ... * S_{r-1}^{a_{r-1}})
-- for all weak compositions a of r-1 into r parts.
-- Returns a HashTable mapping exponent lists to degrees.
mixedIntersections = M -> (
    r := rank M;
    (A, flatToVar) := makeChowRing M;
    al := alphaClass(M, flatToVar);
    Sk := for k from 1 to r-1 list
        sum(select(keys flatToVar, F -> rank(M,F) == r-k), F -> flatToVar#F);
    vals := prepend(al, Sk);
    degOne := al^(r-1);
    R2 := ambient A;
    dc := leadCoefficient lift(degOne, R2);
    deg := f -> (leadCoefficient lift(f + degOne/2, R2)) / dc - 1/2;
    hashTable for e in weakCompositions(r, r-1) list
        e => deg product(r, i -> (vals#i)^(e#i)))

-- Multinomial coefficient n!/(a_1! ... a_k!) for a list of nonnegative integers.
multinomialCoeff = exps -> (
    (sum(exps, i -> i))! / product(exps, i -> i!))

-- Unsigned reduced characteristic coefficients of the interval matroid M[F,G] = (M|G)/F.
-- bar_chi_N(t) = sum_{q=0}^{m-1} (-1)^q mu_q t^{m-1-q}; returns {mu_0, ..., mu_{m-1}}.
intervalMuCoeffs = (M, F, G) -> (
    N := minor(M, F, M.groundSet - G);
    m := rank N;
    if m == 0 then return {1};
    chi := characteristicPolynomial N;
    R := ring chi;
    t := R_0;
    reduced := chi // (t - 1);
    apply(m, q -> (-1)^q * substitute(coefficient(t^(m-1-q), reduced), ZZ)))

-- Closed chain formula for deg(alpha^{a_0} * S_1^{a_1} * ... * S_{r-1}^{a_{r-1}}).
-- This combines Cheng's Chern class formula with the Dastidar-Ross / Eur flag formula
-- for degrees of monomials in the Chow ring. The distinguished element e only enters
-- through alpha; different choices of e give the same class alpha.
-- The inner function takes precomputed data to avoid redundant work.
closedMixedIntersectionInner = (r, d, e, E, flatsByRank, muCache, M, exps) -> (
    total := 0;
    for c in weakCompositions(d, exps#0) do (
        supportRanks := select(toList(1..d), q -> exps#(r - q) + c#(q - 1) > 0);
        if #supportRanks > 0 then (
            dList := apply(supportRanks, q -> exps#(r - q) + c#(q - 1));
            chains := {{} };
            for q in supportRanks do (
                nextChains := {};
                scan(chains, chain -> (
                    prevFlat := if #chain == 0 then {} else chain#(#chain - 1);
                    cands := select(flatsByRank#q, F ->
                        (if #chain == 0 then true else isSubset(prevFlat, F)) and
                        (if c#(q - 1) == 0 then true else member(e, F)));
                    nextChains = join(nextChains, apply(cands, F -> join(chain, {F})))));
                chains = nextChains);
            if #chains > 0 then (
                coeff := multinomialCoeff c;
                chainSum := 0;
                scan(chains, chain -> (
                    k := #chain;
                    partials := {};
                    running := 0;
                    for dj in dList do (
                        running = running + dj;
                        partials = append(partials, running));
                    term := (-1)^(d - k);
                    for j from 0 to k - 1 do (
                        qj := supportRanks#j;
                        dj := dList#j;
                        Dj := partials#j;
                        idx := Dj - qj;
                        nextFlat := if j + 1 < k then chain#(j + 1) else E;
                        key := (chain#j, nextFlat);
                        muList := if muCache#?key then muCache#key
                            else (muCache#key = intervalMuCoeffs(M, chain#j, nextFlat));
                        if idx < 0 or idx > dj - 1 or idx >= #muList then term = 0
                        else term = term * binomial(dj - 1, idx) * muList#idx);
                    chainSum = chainSum + term));
                total = total + coeff * chainSum);
        ));
    total)

closedMixedIntersection = (M, exps) -> (
    r := rank M;
    d := r - 1;
    e := distinguishedElement M;
    E := sort toList(M.groundSet);
    properFlats := select(flats M, F -> #F > 0 and rank(M, F) < r);
    flatsByRank := hashTable for q from 1 to d list
        q => select(properFlats, F -> rank(M, F) == q);
    muCache := new MutableHashTable;
    closedMixedIntersectionInner(r, d, e, E, flatsByRank, muCache, M, exps))

-- Closed chain formula for all mixed intersections appearing in the Chern numbers.
-- Precomputes shared data and memoizes interval characteristic coefficients.
closedMixedIntersections = M -> (
    r := rank M;
    d := r - 1;
    e := distinguishedElement M;
    E := sort toList(M.groundSet);
    properFlats := select(flats M, F -> #F > 0 and rank(M, F) < r);
    flatsByRank := hashTable for q from 1 to d list
        q => select(properFlats, F -> rank(M, F) == q);
    muCache := new MutableHashTable;
    hashTable for exps in weakCompositions(r, r - 1) list
        exps => closedMixedIntersectionInner(r, d, e, E, flatsByRank, muCache, M, exps))

-- Chern numbers from the closed chain formula.
closedChernNumbers = M -> (
    r := rank M;
    if r < 2 then return {};
    (R, ck, formulas) := chernNumberFormula r;
    mv := closedMixedIntersections M;
    b := basis(r - 1, R);
    apply(formulas, (p, poly) -> (
        cs := last coefficients(poly, Monomials => b);
        cn := sum(numcols b, j -> (
            e := flatten exponents b_(0, j);
            lift(cs_(j, 0), QQ) * mv#e));
        (p, cn))))

-- Compare ring-based mixed intersections with the closed chain formula.
verifyClosedMixedIntersections = M -> (
    r := rank M;
    direct := mixedIntersections M;
    closed := closedMixedIntersections M;
    apply(weakCompositions(r, r - 1), e -> (e, direct#e, closed#e, direct#e == closed#e)))

-- Compare ring-based Chern numbers with the closed chain formula.
verifyClosedChernFormula = M -> (
    direct := matroidChernNumbers M;
    closed := closedChernNumbers M;
    apply(#direct, i -> (direct#i#0, direct#i#1, closed#i#1, direct#i#1 == closed#i#1)))

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

-- toddGenus = M -> (
--     r := rank M; n := r - 1;
--     if n == 0 then return 1;
--     (A, flatToVar) := makeChowRing M;
--     al := alphaClass(M, flatToVar);
--     Sk := for k from 1 to r-1 list
--         sum(select(keys flatToVar, F -> rank(M,F) == r-k), F -> flatToVar#F);
--     qc := {1, 1/2, 1/12, 0, -1/720, 0, 1/30240, 0, -1/1209600};
--     Q := x -> sum(min(n+1, #qc), k -> qc#k * x^k);
--     toddClass := product(r-1, k -> Q(Sk#k)) *
--         product(r, i -> Q(al - sum(min(i,r-1), j -> Sk#j)));
--     degOne := al^n; R2 := ambient A;
--     dc := leadCoefficient lift(degOne, R2);
--     deg := f -> (leadCoefficient lift(f + degOne/2, R2)) / dc - 1/2;
--     b := basis(n, A);
--     cs := last coefficients(toddClass, Monomials => b);
--     deg sum(numcols b, j -> (cs_(j,0)) * b_(0,j)))

-- Hirzebruch chi_y genus. For matroid wonderful varieties (rational),
-- h^{p,q} = 0 for p != q, so chi_y = sum (-1)^p * dim(A^p(M)) * y^p.
-- Specializes to: chi_0 = Todd genus = 1, chi_{-1} = Euler char = c_{r-1}.
chiYGenus = M -> (
    r := rank M; n := r - 1;
    (A, flatToVar) := makeChowRing M;
    R := QQ[y];
    sum(n + 1, p -> (-1)^p * (numcols basis(p, A)) * y^p))

-- Kuwata's formula for permutohedral Chern numbers (arXiv:2510.21528 Theorem 1.1).
-- <c_k c_{n-k}, [X_{A_n}]> = (n+1)! * mu_k(n) where
-- mu_k(n) = sum_{j=0}^{floor(k/2)} (1/12)^j * C(k-j,j) * C(n-k-j,j).
kuwataCoeff = (k, n) ->
    (n+1)! * sum(k // 2 + 1, j -> (1/12)^j * binomial(k-j, j) * binomial(n-k-j, j))

-- Boolean specialization of the closed chain formula for I_a(B_r).
-- For the Boolean matroid B_r = U(r,r), every interval [F,G] is again Boolean,
-- so mu_s(B_m) = C(m-1, s). The chain count N_c(q_1,...,q_k) is:
--   if a_0 = 0: r! / (q_1! (q_2-q_1)! ... (r-q_k)!)
--   if a_0 > 0: (t(c)/r) * r! / (q_1! (q_2-q_1)! ... (r-q_k)!)
-- where t(c) = min{q : c_q > 0}.
booleanMixedIntersection = (r, exps) -> (
    d := r - 1;
    a0 := exps#0;
    total := 0;
    for c in weakCompositions(d, a0) do (
        mVals := apply(d, q -> exps#(r - q - 1) + c#q);
        supportRanks := select(toList(1..d), q -> mVals#(q-1) > 0);
        k := #supportRanks;
        if k == 0 then continue;
        dList := apply(supportRanks, q -> mVals#(q-1));
        partials := {};
        running := 0;
        for dj in dList do (running = running + dj; partials = append(partials, running));
        gaps := {};
        prev := 0;
        for q in supportRanks do (gaps = append(gaps, q - prev); prev = q);
        gaps = append(gaps, r - prev);
        chainCount := r! / product(gaps, g -> g!);
        if a0 > 0 then (
            tc := first select(toList(1..d), q -> c#(q-1) > 0);
            chainCount = chainCount * tc / r);
        prod := 1;
        for j from 0 to k - 1 do (
            qj := supportRanks#j;
            mj := mVals#(qj - 1);
            Dj := partials#j;
            idx := Dj - qj;
            qjp1 := if j + 1 < k then supportRanks#(j+1) else r;
            gapSize := qjp1 - qj;
            if idx < 0 or idx > mj - 1 or idx > gapSize - 1 then (prod = 0; break);
            prod = prod * binomial(mj - 1, idx) * binomial(gapSize - 1, idx));
        total = total + multinomialCoeff(c) * chainCount * (-1)^(d - k) * prod;
    );
    total)

booleanMixedIntersections = r -> (
    hashTable for exps in weakCompositions(r, r - 1) list
        exps => booleanMixedIntersection(r, exps))

-- Chern numbers of Boolean matroid (= permutohedral variety) via closed formula.
booleanChernNumbers = r -> (
    if r < 2 then return {};
    (R, ck, formulas) := chernNumberFormula r;
    mv := booleanMixedIntersections r;
    b := basis(r - 1, R);
    apply(formulas, (p, poly) -> (
        cs := last coefficients(poly, Monomials => b);
        cn := sum(numcols b, j -> (
            e := flatten exponents b_(0, j);
            lift(cs_(j, 0), QQ) * mv#e));
        (p, cn))))

-- Hook Chern number: <c_{n-k} c_1^k, [Perm_n]> = (n+1)! * C_{k+1}/(k+1)!
-- where C_m is the m-th Catalan number. Ratio depends only on k, not n.
hookChernNumber = (n, k) ->
    (n+1)! * binomial(2*k+2, k+1) / (k+2)!

-- Klyachko algebra Kly_n = QQ[varpi_1,...,varpi_n] / (varpi_i * alpha_i = 0)
-- where alpha_i = -varpi_{i-1} + 2*varpi_i - varpi_{i+1} (Cartan matrix of A_n).
-- Squarefree monomials form a basis; deg(varpi_1...varpi_n) = 1.
-- Chern numbers: c_lambda = deg(e_{lambda_1}(alpha) * ... * e_{lambda_k}(alpha)).
makeKlyachkoAlgebra = n -> (
    w := getSymbol "w";
    R := QQ[w_1..w_n];
    alpha := apply(n, i -> (
        -1 * (if i > 0 then R_(i-1) else 0_R) +
        2 * R_i +
        -1 * (if i < n-1 then R_(i+1) else 0_R)));
    rels := apply(n, i -> R_i * alpha#i);
    A := R / ideal rels;
    alphaA := apply(alpha, a -> sub(a, A));
    (A, alphaA))

klyachkoChernNumbers = n -> (
    if n < 1 then return {};
    (A, alpha) := makeKlyachkoAlgebra n;
    E := new MutableList from apply(n + 1, k -> if k == 0 then 1_A else 0_A);
    scan(alpha, a -> scan(reverse toList(1..n), k -> E#k = E#k + a * E#(k-1)));
    ck := toList E;
    -- deg(varpi_1...varpi_n) = n! on the permutohedral variety
    degOne := product(n, i -> A_i);
    R2 := ambient A;
    dc := leadCoefficient lift(degOne, R2);
    deg := f -> n! * ((leadCoefficient lift(f + degOne/2, R2)) / dc - 1/2);
    apply(partitions n, p -> (toList p, deg product(toList p, i -> ck#i))))

-- All Chern numbers of the permutohedral variety via equivariant
-- localization (Bott residue formula). For sigma in S_{n+1}, the tangent
-- weights at the fixed point are d_i = sigma(i-1) - sigma(i).
-- c_lambda = sum_{sigma in S_{n+1}} prod_j e_{lambda_j}(d) / prod_i d_i.
permChernNumbers = n -> (
    perms := permutations toList(0..n);
    apply(partitions n, p -> (
        pList := toList p;
        cn := sum(perms, sigma -> (
            d := apply(n, i -> sigma#i - sigma#(i+1));
            es := fold((es, w) -> apply(n+1, j ->
                if j == 0 then es#0 else es#j + es#(j-1) * w),
                prepend(1, toList(n:0)), d);
            product(pList, lam -> es#lam) / last es));
        (pList, cn))))
