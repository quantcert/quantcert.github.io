n := 3;
SympSp := QuantumSymplecticSpace(n);

H, E := Quadrics(SympSp);

quadricsIntersections := { elt[1][3] meet elt[2][3] : elt in CartesianProduct(H,E) };
signatureTypes := [
<7, [ 0,  7, 8 ]>, //1
<7, [ 0,  9, 6 ]>, //2
<6, [ 1,  5, 9 ]>, //3
<5, [ 2,  5, 8 ]>, //4
<5, [ 2,  7, 6 ]>, //5
<4, [ 3,  5, 7 ]>, //6
<3, [ 0,  9, 6 ]>, //7
<3, [ 0, 15, 0 ]>, //8
<3, [ 2,  7, 6 ]>, //9
<3, [ 2,  9, 4 ]>, //10
<3, [ 4,  5, 6 ]>, //11
<3, [ 4,  7, 4 ]>, //12
<3, [ 6,  9, 0 ]>  //13
];
"
The following signatures will be use across the whole script.";
"[<Type,Signature>]:"; 
[ <i,signatureTypes[i]> : i in [1..#signatureTypes] ];

// each of the 1008 element has 15 lines and 15 points :
"
 --- Intersections of quadrics --- ";
"Number of elements:";
#quadricsIntersections;
"Number of lines in each element:";
{ #geometry : geometry in quadricsIntersections };
"Number of points in each element:";
{ #(&join geometry) : geometry in quadricsIntersections };
quadricsSignatures := 
{ <geometry,SHSignature(geometry)> : geometry in quadricsIntersections };
"Signatures as {<{Signature number},Number of occurrences>}:";
{ 
  <
    { i : i in [1..#signatureTypes] | signatureTypes[i] eq signature[2] },
    #[ sig[2] : sig in quadricsSignatures | sig[2] eq signature[2] ]
  > : signature in quadricsSignatures 
};

type9 := { signature[1] : signature in quadricsSignatures | 
  signature[2] eq signatureTypes[9] };
type9A := 
{ 
  geometry 
    : geometry in type9 
    | InnerProduct(
        [ point : point in &join geometry | NumberOfI(point) eq 2][1], 
        [ point : point in &join geometry | NumberOfI(point) eq 2][2]
      ) eq 0
};
// type9B := type9 diff type9A;
type9B := 
{ 
  geometry 
    : geometry in type9 
    | InnerProduct(
        [ point : point in &join geometry | NumberOfI(point) eq 2][1], 
        [ point : point in &join geometry | NumberOfI(point) eq 2][2]
      ) eq 1
};
"Number of geometries with a signature of type 9 such that the operators
containing the identity twice do not commute (called type 9A):";
#type9A;
"Number of geometries with a signature of type 9 such that the operators
containing the identity twice commute (called type 9B):";
#type9B;

// -------------------------------------------------------------------------- //

P := PerpSets(SympSp);

// we generate from the points and not from the lines because ??? TODO
perpsetsPointsIntersections := 
{ elt[1][2] meet elt[2][2] : elt in CartesianProduct(P,P) |
  InnerProduct(elt[1][1],elt[2][1]) eq 1 };
lines := IsotropicSubspaces(SympSp,1);
perpsetsLinesIntersections := 
{ 
  { line : line in lines | line subset geometryPoints } : 
  geometryPoints in perpsetsPointsIntersections 
};

// each of the 336 element has 15 lines and 15 points :
"
 --- Intersections of perpsets --- ";
"Number of elements:";
#perpsetsLinesIntersections;
"Number of lines in each element:";
{ #geometry : geometry in perpsetsLinesIntersections };
"Number of points in each element:";
{ #(&join geometry) : geometry in perpsetsLinesIntersections };
perpsetsSignatures := 
{ 
  <geometry,SHSignature(geometry)> 
    : geometry in perpsetsLinesIntersections 
};
"Signatures as {<{Signature number},Number of occurrences>}:";
{ 
  <
    { i : i in [1..#signatureTypes] | signatureTypes[i] eq signature[2] },
    #[ sig[2] : sig in perpsetsSignatures | sig[2] eq signature[2] ]
  > : signature in perpsetsSignatures 
};
 
// -------------------------------------------------------------------------- //

doilies := perpsetsLinesIntersections join quadricsIntersections;
ovoidsMap := &join
{ { <doily,ovoid> : ovoid in InnerOvoids(doily) } : doily in doilies };
"
 --- Properties of W(5,2)'s doilies --- ";
"Number of doilies containing a given ovoid:";
{ 
  #{ doily : doily in doilies | <doily,ovoid> in ovoidsMap } 
    : ovoid in { ovoid[2] : ovoid in ovoidsMap } 
};

gridsMap := &join
{ { <doily,grid> : grid in InnerGrids(doily) } : doily in doilies };
"Number of doilies containing a given grid:";
{ 
  #{ doily : doily in doilies | <doily,grid> in gridsMap } 
    : grid in { grid[2] : grid in gridsMap } 
};

perpsetsMap := &join
{ { <doily,perpset> : perpset in InnerPerpsets(doily) } : doily in doilies };
"Number of doilies containing a given perpset:";
{ 
  #{ doily : doily in doilies | <doily,perpset> in perpsetsMap } 
    : perpset in { perpset[2] : perpset in perpsetsMap } 
};

// -------------------------------------------------------------------------- //

allSignatures := quadricsSignatures join perpsetsSignatures;
SympSp2 := QuantumSymplecticSpace(2);
Doily2 := WLines(SympSp2)[1][3];
ovoids := Elliptics(SympSp2);
increasedWithOvoidsDoilies := 
&join
{ 
  &join
  { 
    {
      IncreaseGeometry(Doily2,{<pos,operator,ovoid[2]>}) : pos in [1..n]
    } : ovoid in ovoids
  } : operator in Subsequences({0,1},2) | operator ne [0,0]
};
hyperplans := Elliptics(SympSp2) join PerpSets(SympSp2) join Hyperbolics(SympSp2);
increasedWithHyperplansDoilies := 
&join
{ 
  &join
  { 
    {
      IncreaseGeometry(Doily2,{<pos,operator,ovoid[2]>}) : pos in [1..n]
    } : ovoid in hyperplans
  } : operator in Subsequences({0,1},2) | operator ne [0,0]
};
"
 --- Doilies of W(5,2) generated from W(3,2) --- ";
"All doilies with a signature of type 1 to 7 are genuine W(5,2) doilies:";
not 
{ 
  signature[1] : signature in allSignatures | signature[2] in signatureTypes[1..7] 
} subset increasedWithHyperplansDoilies;
"27 of the 36 doilies with a signature of type 8 are genuine W(5,2) doilies:";
#
{ 
  signature[1] 
    : signature in allSignatures 
    | (signature[2] eq signatureTypes[8]) 
        and 
      (not signature[1] in increasedWithHyperplansDoilies)
}
eq 27;
"All doilies with a signature of type 9A are genuine W(5,2) doilies:";
IsEmpty(increasedWithHyperplansDoilies meet type9A);
"All doilies with a signature of type 9B are extended W(3,2)'s:";
type9B eq increasedWithOvoidsDoilies;
"All doilies with a signature of type 10 to type 12 are extended W(3,2)'s:";
{ 
  signature[1] 
    : signature in allSignatures 
    | signature[2] in signatureTypes[10..12] 
} 
subset increasedWithHyperplansDoilies;
"All doilies with a signature of type 13 are extended W(3,2)'s:";
&and
{ 
  &and{ 
    signature[2] eq signatureTypes[13] 
      : signature in allSignatures 
      | signature[1] eq IncreaseGeometry(Doily2,{<pos,[0,0],{}>})
  } : pos in [1..n] 
};

// -------------------------------------------------------------------------- //

"
 --- Deep and zero-points of W(5,2) --- ";

type3 := { signature[1] : signature in allSignatures | 
  signature[2] eq signatureTypes[3] };

"The deep and zero-points of doilies of type 3 of W(5,2) are forming a
tricentric triad, i.e. are isomorphic to W(1,2)";
"Number of points in those tricentric triads:";
{ #(DeepPoints(doily) join ZeroPoints(doily)) : doily in type3 };

// -------------------------------------------------------------------------- //

"
 --- Finer study of contextuality --- ";
"Contextuality degree per signature type";
allSignaturesDegrees :=
{
  < signature[1], signature[2], ContextualityDegree(signature[1]) > 
    : signature in allSignatures 
};
sigDegGeometries :=
[
  < 
    sigDeg, 
    { 
      geomSigDeg[1] 
        : geomSigDeg in allSignaturesDegrees
        | geomSigDeg[2] eq signatureTypes[sigDeg[1]] 
            and 
          geomSigDeg[3] eq sigDeg[2]
    }
  > : sigDeg in CartesianProduct([1..#signatureTypes],[1..3])
];
"[<<Signature type,degree>, Number of such combinations, Example>]:
(The examples may be disabled)";
[ 
  < 
    sigDegGeometry[1], 
    #sigDegGeometry[2]
    // , not #sigDegGeometry[2] eq 0 select Random(sigDegGeometry[2]) else {}
  > 
    : sigDegGeometry in sigDegGeometries 
];

// -------------------------------------------------------------------------- //

"
 --- Study of points order per signature --- ";
"Order of the points of the geometries of each signature:";
linesValues := { <line,SetSign(line,PauliOperatorCanonical)> : line in Lines(SympSp) };
[
  < 
    sigIndex, 
    {
      PointsOrders(geomSig[1],linesValues)
        : geomSig in allSignatures
        | geomSig[2] eq signatureTypes[sigIndex]
    } 
  >
    : sigIndex in [1..#signatureTypes]
];

// exit;
