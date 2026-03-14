load "main.m2";

-- Matroid Chern numbers (Cheng, arXiv:2510.06609)
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

-- Strict transform Chern numbers at each blow-up stage, using Sect. 4 of arXiv:2510.06609
strictTransformChernTable specificMatroid "fano"
strictTransformChernTable uniformMatroid(3,5)
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

-- Kuwata's formula for permutohedral Chern numbers (arXiv:2510.21528)
-- <c_k c_{n-k}, [X_{A_n}]> = (n+1)! * mu_k(n)
kuwataCoeff(1, 2)  -- c_1^2 on X_(A_2) = 6
kuwataCoeff(1, 3)  -- c_1*c_2 on X_(A_3) = 24
kuwataCoeff(2, 4)  -- c_2^2 on X_(A_4) = 130
kuwataCoeff(1, 5)  -- c_1*c_4 on X_(A_5) = 720

-- Closed formula tests for Chern numbers
closedChernNumbers uniformMatroid(3, 3)
closedChernNumbers uniformMatroid(3, 5)
closedChernNumbers(matroid completeGraph 4)
closedChernNumbers(specificMatroid "fano")
closedChernNumbers(specificMatroid "pappus")
closedChernNumbers uniformMatroid(4, 4)
closedChernNumbers uniformMatroid(4, 6)
closedChernNumbers(matroid completeGraph 5)