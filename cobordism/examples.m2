load "main.m2";

-- Matroid Chern numbers (arXiv:2510.06609)
matroidChernNumbers uniformMatroid(3, 3)
matroidChernNumbers uniformMatroid(4, 4)
matroidChernNumbers(matroid completeGraph 4)
matroidChernNumbers(specificMatroid "fano")
matroidChernNumbers uniformMatroid(3, 5)

-- Toric Chern numbers
toricChernNumbers permutohedralVariety 2
toricChernNumbers permutohedralVariety 3
toricChernNumbers toricProjectiveSpace 3

chernTable 3
chernTable 4

-- Strict transform Chern numbers at each blow-up stage
strictTransformChernTable specificMatroid "fano"
strictTransformChernTable uniformMatroid(3,4)
strictTransformChernTable matroid completeGraph 4
strictTransformChernTable matroid completeGraph 5

-- Universal Chern number formula (rank 3 and 4)
(R3, ck3, f3) = chernNumberFormula 3
scan(f3, (p,poly) -> print(toString p | " : " | toString poly))
(R4, ck4, f4) = chernNumberFormula 4
scan(f4, (p,poly) -> print(toString p | " : " | toString poly))

-- Mixed intersection numbers
mixedIntersections uniformMatroid(3, 3)
mixedIntersections(matroid completeGraph 4)

-- Verify combinatorial formula matches direct computation
combinatorialChernNumbers uniformMatroid(3, 3)
combinatorialChernNumbers uniformMatroid(4, 4)
combinatorialChernNumbers(matroid completeGraph 4)
combinatorialChernNumbers(specificMatroid "fano")

-- Todd genus (always 1 for matroid wonderful varieties)
toddGenus uniformMatroid(3, 3)
toddGenus uniformMatroid(4, 4)
toddGenus(matroid completeGraph 4)

-- chi_y genus: polynomial in y with Betti number coefficients
-- chi_0 = 1 (Todd), chi_{-1} = Euler characteristic
chiYGenus uniformMatroid(3, 3)
chiYGenus uniformMatroid(4, 4)
chiYGenus uniformMatroid(5, 5)
chiYGenus(matroid completeGraph 4)
chiYGenus(matroid completeGraph 5)