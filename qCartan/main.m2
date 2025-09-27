loadPackage "Polyhedra";

R = ZZ;

Cartan = (n,q) -> (
M = matrix table(n, n, (i,j) -> if j == i+1 then 1_R else 0_R); 
(q+1)*id_(R^n) - M - q * transpose M
)


partitionCone = method()
partitionCone(List,List,Matrix,Matrix) := (J,K,A,B) -> (
     if (numcols A != numcols B) then error("A and B should have the same number of columns.");
     n = numcols A;
     if (isSubset(set J, set toList(0..n-1)) == false) then error("J should be a subset of 0..n-1");
     if (isSubset(set K, set toList(0..n-1)) == false) then error("K should be a subset of 0..n-1");
     if (#intersect(set J,set K) > 0) then error("J and K should be disjoint");
    coneFromVData (A_J | -(B_K))
)

partitionFan = method()
partitionFan := (n,q) -> (
    if (n <= 0) then error("n should be a positive integer.");
    if (q <= 0) then error("q should be a positive integer.");
    R = ZZ;
    A = Cartan(n,q);
    B = id_(R^n);
    F = {};
    for k from 0 to n do (
        for J in subsets(toList(0..n-1),k) do (
           K = set toList(0..n-1) - set J;
           F = append(F, {J, toList K});
        );
    );
    L = {};
    for pair in F do (
        J = pair#0;
        K = pair#1;
        L = append(L, partitionCone(J,K,A,B));
    );
    fan L
)



SigmaNQ = partitionFan(2,3)

PNQ = polytope SigmaNQ

vertices PNQ
