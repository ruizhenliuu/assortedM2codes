needsPackage "Matroids"
-- Translation of Eulerchar.sage to Macaulay2
-- Ground set should be numbers and should contain 0

--
-- Functions to create classes in A(M)
--

-- Creates the Chow ring of M and a HashTable mapping flats to ring variables
-- Returns (A, flatToVar)
makeChowRing = M -> (
    I := idealChowRing M;
    R := ring I;
    A := R/I;
    properFlats := select(flats M, F -> #F > 0 and rank(M,F) < rank M);
    flatToVar := hashTable apply(#properFlats, i -> properFlats#i => A_i);
    (A, flatToVar)
)

-- Creates the divisor alpha = sum x_F for flats F containing 0
alphaClass = flatToVar -> (
    sum(select(keys flatToVar, F -> member(0, F)), F -> flatToVar#F)
)

-- Computes the Chern roots of S_N restricted to A(M)
chernRootSN = (M, flatToVar, N) -> (
    al := alphaClass flatToVar;
    pf := keys flatToVar;
    for i from 0 to rank(N) - 1 list (
        al - sum(select(pf, F -> rank(N, F) > i), F -> flatToVar#F)
    )
)

-- Returns the first Chern class of the line bundle P(N), in A(M)
matroidToChow = (M, flatToVar, N) -> (
    Nd := dual N;
    al := alphaClass flatToVar;
    c1 := -(rank Nd) * al;
    for F in keys flatToVar do (
        c1 = c1 + (rank(Nd, F)) * flatToVar#F;
    );
    c1
)

-- Computes Chern roots of Q_N restricted to A(M)
chernRootQN = (M, flatToVar, N) -> (
    Nd := dual N;
    nCorank := #(N.groundSet) - rank N;
    trunclasses := for i from 0 to nCorank list (
        c := matroidToChow(M, flatToVar, Nd);
        if i < nCorank then Nd = truncate Nd;
        c
    );
    for i from 0 to nCorank - 1 list (
        trunclasses#i - trunclasses#(i+1)
    )
)

-- Takes a list of Chern roots, returns the Chern classes
chernClasses = (L, A) -> (
    if #L == 0 then return {1_A};
    prepend(1_A, for k from 1 to #L list (
        sum(subsets(L, k), term -> product term)
    ))
)

--
-- Functions related to polynomials
--

-- Eulerian polynomial A_k(t)
eulerianPoly = (k, R) -> (
    t := R_0;
    if k <= 1 then return 1_R;
    a := {1};
    for n from 2 to k do (
        a = for j from 0 to n-1 list (
            (j+1) * (if j < #a then a#j else 0) +
            (n-j) * (if j-1 >= 0 then a#(j-1) else 0)
        );
    );
    sum(k, j -> (a#j) * t^j)
)

-- Lagrange interpolation through points {{x0,y0}, ...}
lagrangeInterp = (pts, R) -> (
    n := #pts;
    V := sub(matrix apply(n, i -> apply(n, j -> (pts#i#0)^j)), QQ);
    Y := sub(matrix apply(n, i -> {pts#i#1}), QQ);
    c := solve(V, Y);
    t := R_0;
    sum(n, i -> (c_(i,0)) * t^i)
)

-- Converts Ehrhart polynomial P(t) to h* vector
ehrhartToHstar = (P, R) -> (
    t := R_0;
    d := (degree P)#0;
    coeffList := apply(d + 1, i -> coefficient(t^i, P));
    ans := (1-t)^d * coeffList#0;
    for i from 1 to d do (
        ans = ans + coeffList#i * t * (1-t)^(d-i) * eulerianPoly(i, R);
    );
    apply(d + 1, i -> coefficient(t^i, ans))
)

-- h-vector to f-vector for a (d-1)-dimensional complex
htofvector = (L, d) -> (
    h := L | apply(d + 1 - #L, i -> 0);
    prepend(1, for i from 1 to d list (
        sum(i + 1, k -> binomial(d - k, i - k) * h#k)
    ))
)

-- f-vector to h-vector
ftohvector = L -> (
    d := #L - 1;
    for k from 0 to d list (
        sum(k + 1, i -> (-1)^(k-i) * binomial(d - i, k - i) * L#i)
    )
)

--
-- Functions to compute Euler characteristics
--

-- Computes chi(M, [P(N)]^num)
-- Uses the HRR-type formula: deg_M((1 + alpha + ...)c(S_{N^*}^vee)^num)
eulerChar = (M, N, num) -> (
    (A, flatToVar) := makeChowRing M;
    R := ambient A;
    Nd := dual N;
    chernSList := chernClasses(chernRootSN(M, flatToVar, Nd), A);
    totalclass := sum(#chernSList, i -> (-1)^i * chernSList#i);
    al := alphaClass flatToVar;
    totalalpha := sum(rank M, i -> al^i);
    degOneElt := al^(rank M - 1);
    p := totalclass^num * totalalpha;
    -- degree extraction via leading coefficient trick
    pplus := p + (1/2) * degOneElt;
    (leadCoefficient lift(pplus, R)) / (leadCoefficient lift(degOneElt, R)) - 1/2
)

-- Compute the h* vector of P(N) in A(M)
hstarMN = (M, N) -> (
    R := QQ[symbol t];
    pts := prepend({0, 1}, for i from 1 to rank M - 1 list {i, eulerChar(M, N, i)});
    ehr := lagrangeInterp(pts, R);
    d := (degree ehr)#0;
    ans := ehrhartToHstar(ehr, R);
    ans | apply(d + 1 - #ans, i -> 0)
)

-- Compute the Ehrhart polynomial of P(N) on M
ehrhartPoly' = (M, N) -> (
    R := QQ[symbol t];
    pts := prepend({0, 1}, for i from 1 to rank M - 1 list {i, eulerChar(M, N, i)});
    lagrangeInterp(pts, R)
)

-- Computes chi(M, P(N) tensor P(P)^{-1})
eulerMatroidPolytopeDiff = (M, N, P) -> (
    (A, flatToVar) := makeChowRing M;
    R := ambient A;
    al := alphaClass flatToVar;
    totalalpha := sum(rank M, i -> al^i);
    degOneElt := al^(rank M - 1);
    Nd := dual N;
    Pd := dual P;
    chernSList := chernClasses(chernRootSN(M, flatToVar, Nd), A);
    chernQList := chernClasses(chernRootQN(M, flatToVar, Pd), A);
    totalclassS := sum(#chernSList, i -> (-1)^i * chernSList#i);
    totalclassQ := sum(#chernQList, i -> (-1)^i * chernQList#i);
    p := totalclassS * totalclassQ * totalalpha;
    pplus := p + (1/2) * degOneElt;
    (leadCoefficient lift(pplus, R)) / (leadCoefficient lift(degOneElt, R)) - 1/2
)



