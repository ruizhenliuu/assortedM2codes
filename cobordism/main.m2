needsPackage "NormalToricVarieties"

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

permutohedralVariety = n -> (
    e := append(apply(n, i -> apply(n, j -> if i == j then 1 else 0)), apply(n, j -> -1));
    (rys, subs) := (apply(delete({}, delete(toList(0..n), subsets toList(0..n))),
        S -> sum(S, i -> e#i)),
        delete({}, delete(toList(0..n), subsets toList(0..n))));
    idx := hashTable apply(#rys, i -> rys#i => i);
    normalToricVariety(rys, unique apply(permutations(n + 1),
        sigma -> sort apply(n, k -> idx#(sum(sort take(sigma, k + 1), i -> e#i))))))
