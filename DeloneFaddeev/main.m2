-- working with ZZ/32003 to get dimension,
R = ZZ/32003[
  c111, c112, c113, c121, c122, c123, c131, c132, c133,
  c211, c212, c213, c221, c222, c223, c231, c232, c233,
  c311, c312, c313, c321, c322, c323, c331, c332, c333,
  d1, d2, d3
];

use R;

c = (i,j,k) -> value getSymbol concatenate {"c", toString i, toString j, toString k};

commute = for i from 1 to 3 list for j from i+1 to 3 list (
  for k from 1 to 3 list (c(i,j,k) - c(j,i,k))
);

Icomm = ideal flatten flatten commute;


assoc = for i from 1 to 3 list for j from 1 to 3 list for l from 1 to 3 list (
  for m from 1 to 3 list (
    (sum for k from 1 to 3 list (c(i,j,k)*c(k,l,m))) -
    (sum for k from 1 to 3 list (c(j,l,k)*c(i,k,m)))
  )
);

Iassoc = ideal flatten flatten flatten assoc;

#(flatten flatten flatten assoc) -- 81


d = k -> value getSymbol concatenate {"d", toString k};

unit = for i from 1 to 3 list for m from 1 to 3 list (
  sum(for j from 1 to 3 list (d(j)*c(j,i,m))) - (if i==m then 1 else 0)
);
Iunit = ideal flatten flatten unit;

I = Icomm + Iassoc + Iunit

mingens I 

dim (R/I) -- 9

Icub = I + ideal(d1-1,d2,d3)

dim (R/Icub) -- 6

InormCub = Icub + ideal(c223,c323)

dim (R/InormCub) -- 4 

-- Use Groebner basis to find a parameterisation 

gens gb InormCub

-*
The output of gens gb InormCub is 

o20 = | d3 d2 d1-1 c323 c321 c313-1 c312 c311 c233 c232-c322 c231 c223 c221
      -------------------------------------------------------------------------
      c213 c212-1 c211 c133-1 c132 c131 c123 c122-1 c121 c113 c112 c111-1
      -------------------------------------------------------------------------
      c322^2-c222c332-c322c333-c331 |
*-

gens gb Icub

-*
The output of gens gb Icub is
 | d3 d2 d1-1 c313-1 c312 c311 c233-c323 c232-c322 c231-c321 c213 c212-1 c211 c133-1 c132 c131
      -----------------------------------------------------------------------------------------------
      c123 c122-1 c121 c113 c112 c111-1 c322c323-c223c332+c321 c322^2-c222c332+c323c332-c322c333-c331
      -----------------------------------------------------------------------------------------------
      c321c322+c323c331-c221c332-c321c333 c223c322-c222c323+c323^2-c223c333-c221
      -----------------------------------------------------------------------------------------------
      c222c321-c221c322-c321c323+c223c331 c323^2c331+c223c321c332-c221c323c332-c321c323c333-c321^2
      -----------------------------------------------------------------------------------------------
      c222c323^2-c323^3-c223^2c332+c223c323c333+c223c321+c221c323 |
*-
