W3 := QSS(3);

// pair : 2i -> n+i ; j -> j/2+n
// impair : 2i-1 -> i ; j -> (j+1)/2
function basisChange(i,n)
  return Integers()!(IsEven(i) select (i+n)/2 else (i+1)/2);
end function;

function pointBasisChange(p)
  n := #Eltseq(p);
  S := Parent(p);
  return S![ p[basisChange(i,n)] : i in [1..n] ];
end function;

// lt n : i -> 2i-1
// gt n : i -> 2(i-n)
function invBasisChange(i,n)
  return Integers()!(i lt n/2+1 select 2*i-1 else 2*(i-n/2));
end function;

function pointInvBasisChange(p)
  n := #Eltseq(p);
  S := Parent(p);
  return S![ p[invBasisChange(i,n)] : i in [1..n] ];
end function;

function seventhCoord(p) 
// p should be in the basis of the article, or in the long time we should 
// convert all coordinates to ours
  return p[1]*p[4]+p[2]*p[5]+p[3]*p[6];
end function;

function Q(p)
  p1 := pointBasisChange(p);
  p_7 := seventhCoord(p1);
  return p1[1]*p1[4]+p1[2]*p1[5]+p1[3]*p1[6]+p_7;
end function;

cayleyHexagonPoints := { p : p in W3 | IsZero(Q(p)) and not IsZero(p) };

function com(x,y,i,j)
  return x[i]*y[j]-x[j]*y[i];
end function;

// TODO, deduce this from Q
function pointsAligned(x,y)
  xp := pointInvBasisChange(x);
  yp := pointInvBasisChange(y);
  xe := [ c : c in Eltseq(xp) ] cat [ seventhCoord(xp) ];
  ye := [ c : c in Eltseq(yp) ] cat [ seventhCoord(yp) ];
  return 
    com(xe,ye,6,2) eq com(xe,ye,1,7) and 
    com(xe,ye,1,3) eq com(xe,ye,7,2) and 
    com(xe,ye,2,4) eq com(xe,ye,3,7) and 
    com(xe,ye,3,5) eq com(xe,ye,7,4) and 
    com(xe,ye,4,6) eq com(xe,ye,5,7) and 
    com(xe,ye,5,1) eq com(xe,ye,7,6) and 
    com(xe,ye,1,4)+com(xe,ye,2,5)+com(xe,ye,3,6) eq 0;
end function;

cayleyHexagon := 
{ 
  l 
    : l in Lines(W3) 
    | &and{ pointsAligned(xy[1],xy[2]) : xy in Subsequences(l,2) } 
};

function epsilon(p)
  S := Parent(p);
  p1 := pointInvBasisChange(p);
  p_7 := seventhCoord(p1);
  p2 := 
  [ 
    p1[1] + p1[6] + p1[4]*p1[6] + p_7*p1[5],
    p1[2] + p1[3] + p1[3]*p1[5] + p_7*p1[4],
    p1[3], p1[4], p1[5], p1[6]
  ];
  return S!pointBasisChange(p2);
end function;

skewCayleyHexagon := { { epsilon(p) : p in l } : l in cayleyHexagon };

for el in Elliptics(W3) do
  el[1];
  newGeom := Lines(W3) diff (el[3] meet skewCayleyHexagon);
  #newGeom; //306
  IsContextual(newGeom); //true
end for;

for hyp in Hyperbolics(W3) do
  hyp[1];
  newGeom := Lines(W3) diff (hyp[3] meet skewCayleyHexagon);
  #newGeom; //294
  IsContextual(newGeom); //true
end for;

for p in { p : p in W3 | not IsZero(p) } do
  p;
  pts, ls := Quadric(p);
  newGeom := Lines(W3) diff (ls meet cayleyHexagon);
  #newGeom;
  IsContextual(newGeom);
end for;

"num of lines where sum is zero";
#{ l : l in skewCayleyHexagon | &+l eq W3!0 };
"num of lines where points in line commute";
#{ 
  l 
    : l in skewCayleyHexagon 
    | &and{ IsZero(InnerProduct(pair[1],pair[2])) : pair in Subsequences(l,2) } 
};
// TODO: check contextuality degree of the complement of this skew caylex hex

W3Lines := Lines(W3);
oldCard := #skewCayleyHexagon;
copy := skewCayleyHexagon;
for i in [1..5] do
  while #copy lt oldCard+i do
    copy join:= { Random(W3Lines) };
  end while;
  printf "skew Caley hexagon with %o more random line(s) is contextual ?\n", i;
  IsContextual(copy);
end for;

function transvection(p)
  function Tp(q)
    return q+InnerProduct(p,q)*p;
  end function;
  return Tp;
end function;

// identityPoints := {};
transformationClassic := {};
for points in CartesianPower(W3,2) do
  transvections := [ transvection(p) : p in points ];
  classicCopy := cayleyHexagon;
  for T in transvections do
    classicCopy := { { T(x) : x in l } : l in classicCopy };
  end for;
  transformationClassic join:= { classicCopy };
end for;

transformationSkew := {};
for points in CartesianPower(W3,3) do
  transvections := [ transvection(p) : p in points ];
  skewCopy := skewCayleyHexagon;
  for T in transvections do
    skewCopy := { { T(x) : x in l } : l in skewCopy };
  end for;
  transformationSkew join:= { skewCopy };
  // if skewCopy eq skewCayleyHexagon then
  //   identityPoints join:= { points };
  // end if;
end for;

"Number of skew embedding obtained:";
#transformationSkew;
"Number of classical embedding obtained:";
#transformationClassic;
"Are the complement of each one of them contextual ?";
{ IsContextual(W3Lines diff hexa) : hexa in transformationSkew };
{ IsContextual(W3Lines diff hexa) : hexa in transformationClassic };

p := Random({ e : e in W3 | not IsZero(e) });

cayleyHyperplane := PerpSet(p) meet cayleyHexagon;
cayleyHyperplaneComplement1 := W3Lines diff cayleyHyperplane;
cayleyHyperplaneComplement1 := W3Lines;

skewCayleyHyperplane := PerpSet(p) meet skewCayleyHexagon;
skewCayleyHyperplaneComplement1 := W3Lines diff skewCayleyHyperplane;

// le choix de epsilon dépend de la quadrique qu'on prend, ici la transfo ne
// marche pas car on n'a pas la bonne quadrique.

P := PerpSets(W3);

perpsetsPointsIntersections := 
{ elt[1][2] meet elt[2][2] : elt in CartesianProduct(P,P) |
  InnerProduct(elt[1][1],elt[2][1]) eq 1 };
lines := QuantumSubspaces(SympSp,1);
perpsetsLinesIntersections := 
{ 
  { line : line in lines | line subset geometryPoints } : 
  geometryPoints in perpsetsPointsIntersections 
};

H, E := Quadrics(SympSp);

quadricsIntersections := { elt[1][3] meet elt[2][3] : elt in CartesianProduct(H,E) };

doilies := perpsetsLinesIntersections join quadricsIntersections;

chexCompl := Lines(W3) diff caleyHexagon;
"How many doilies in the complement of the classis chex ?";
#{ doily : doily in doilies | doily subset chexCompl };

SchexCompl := Lines(W3) diff skewCayleyHexagon;
"How many doilies in the complement of the skew chex ?";
#{ doily : doily in doilies | doily subset SchexCompl };

// TODO : if no results from the lines above, do the same work for the mermin 
// squares

// mSquare := ;

// "How many Mermin squares in the complement of the classis chex ?";
// #{ mSquare : mSquare in mSquares | mSquare subset chexCompl };

// "How many Mermin squares in the complement of the skew chex ?";
// #{ mSquare : mSquare in mSquares | mSquare subset SchexCompl };
