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

function transvection(p)
  function Tp(q)
    return q+InnerProduct(p,q)*p;
  end function;
  return Tp;
end function;

transformationClassic := {};
transformationSkew := {};
// identityPoints := {};
for points in CartesianPower(W3,4) do
  transvections := [ transvection(p) : p in points ];
  skewCopy := skewCayleyHexagon;
  // classicCopy := cayleyHexagon;
  for T in transvections do
    skewCopy := { { T(x) : x in l } : l in skewCopy };
    // classicCopy := { { T(x) : x in l } : l in classicCopy };
  end for;
  // transformationClassic join:= { classicCopy };
  transformationSkew join:= { skewCopy };
  // if skewCopy eq skewCayleyHexagon then
  //   identityPoints join:= { points };
  // end if;
end for;

"Number of skew embedding obtained:";
#transformationSkew;
// #transformationClassic;
"Are the complement of each one of them contextual ?";
{ IsContextual(W3Lines diff hexa) : hexa in transformationSkew };
// { IsContextual(W3Lines diff hexa) : hexa in transformationClassic };
