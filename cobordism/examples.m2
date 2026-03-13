load "main.m2"

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
