signatureTypes := [ <130,[3,9,33,18]>,
                    <126,[0,24,0,39]>,
                    <126,[1,13,27,22]>,
                    <126,[2,10,30,21]>,
                    <122,[1,15,27,20]>,
                    <122,[2,10,30,21]>,
                    <118,[0,16,32,15]>,
                    <118,[3,9,33,18]>,
                    <118,[3,11,25,24]>,
                    <114,[1,15,27,20]>,
                    <114,[1,17,27,18]>,
                    <114,[3,13,25,22]>,
                    <114,[4,12,28,19]>,
                    <110,[3,15,25,20]>,
                    <110,[5,11,23,24]>,
                    <106,[5,13,23,22]>,
                    <102,[1,21,27,14]>,
                    <102,[2,18,30,13]>,
                    <102,[3,15,25,20]>,
                    <102,[4,12,28,19]>,
                    <90,[0,36,0,27]>,
                    <90,[2,22,30,9]>,
                    <90,[3,9,33,18]>,
                    <90,[3,21,25,14]>,
                    <90,[4,16,28,15]>,
                    <90,[5,15,31,12]>,
                    <90,[6,18,26,13]>,
                    <90,[7,17,21,18]>,
                    <90,[9,27,27,0]> ];

n := 4;
SympSp := QuantumSymplecticSpace(n);

H, E := Quadrics(SympSp);

quadricsIntersections := { elt[1][3] meet elt[2][3] : elt in CartesianProduct(H,E) };

"
The following signatures will be use across the whole script.";
"[<Type,Signature>]:"; 
[ <i,signatureTypes[i]> : i in [1..#signatureTypes] ];

"
 --- Intersections of quadrics --- ";
"Number of elements:";
#quadricsIntersections;
"Number of lines in each element:";
{ #geometry : geometry in quadricsIntersections };
"Number of points in each element:";
{ #(&join geometry) : geometry in quadricsIntersections };
"Are they all contextual ?";
&and{ IsContextual(geometry) : geometry in quadricsIntersections };
"Signatures:";
quadricsSignatures := [ SHSignature(geometry) : geometry in quadricsIntersections ];
{ <signature,#[ sig : sig in quadricsSignatures | sig eq signature ]> : 
  signature in Seqset(quadricsSignatures) };

// -------------------------------------------------------------------------- //

P := PerpSets(SympSp);

// we generate from the points and not from the lines because ??? TODO
perpsetsPointsIntersections := 
{ elt[1][2] meet elt[2][2] : elt in CartesianProduct(P,P) |
  InnerProduct(elt[1][1],elt[2][1]) eq 1 };
lines := QuantumSubspaces(SympSp,1);
perpsetsLinesIntersections := 
{ 
  { line : line in lines | line subset geometryPoints } : 
  geometryPoints in perpsetsPointsIntersections 
};

"
 --- Intersections of perpsets --- ";
"Number of elements:";
#perpsetsLinesIntersections;
"Number of lines in each element:";
{ #geometry : geometry in perpsetsLinesIntersections };
"Number of points in each element:";
{ #(&join geometry) : geometry in perpsetsLinesIntersections };
"Are they all contextual ?";
&and{ IsContextual(geometry) : geometry in quadricsIntersections };
"Signatures:";
perpsetsSignatures := [ SHSignature(geometry) : geometry in perpsetsLinesIntersections ];
{ <signature,#[ sig : sig in perpsetsSignatures | sig eq signature ]> : 
  signature in Seqset(perpsetsSignatures) };

// -------------------------------------------------------------------------- //
"
 --- Deep points and zero-points --- ";
w52s := perpsetsLinesIntersections join quadricsIntersections;
"Number of doilies in W(7,2)";
#w52s;
doiliesPoints := 
{ 
  DeepPoints(geometry) join ZeroPoints(geometry) : geometry in w52s 
};
"Number of non empty sets of deep points and zero-points amongst those doilies:";
#{ pointSet : pointSet in doiliesPoints | not IsZero(#pointSet) };
"Number of doilies composed of deep points and zero-points amongst those doilies:";
#{ pointSet : pointSet in doiliesPoints | #pointSet eq 15 };
"Cardinal of sets of deep and zero-points:";
{ #pointSet : pointSet in doiliesPoints };

// -------------------------------------------------------------------------- //
"
 --- Properties of W(7,2)'s doilies --- ";
SympSp2 := QuantumSymplecticSpace(2);
Doily2 := WLines(SympSp2)[1][3];
hyperplanesPoints := 
{
  geometry[2] 
    : geometry in 
      Elliptics(SympSp2) join 
      Hyperbolics(SympSp2) join 
      PerpSets(SympSp2)
};
"Number of hyperplanes in W(3,2):";
#hyperplanesPoints;
PauliOperators := { operator : operator in Subsequences({0,1},2) | operator ne [0,0] };
posOpNeutrals1 := CartesianProduct(<[1..n],PauliOperators,hyperplanesPoints>);
posOpNeutrals2 := 
{ 
  posOpNeutral 
    : posOpNeutral in CartesianPower(posOpNeutrals1,2)
    | not posOpNeutral[1][1] eq posOpNeutral[2][1]
};
increased2Doilies := 
{
  IncreaseGeometry(Doily2,{posOpNeutral[1],posOpNeutral[2]})
    : posOpNeutral in posOpNeutrals2 
};
"Number of doilies of W(7,2) resulting from those hyperplanes:";
#increased2Doilies;

SympSp3 := QuantumSymplecticSpace(3);
H, E := Quadrics(SympSp3);
quadrics3Intersections := { elt[1][3] meet elt[2][3] : elt in CartesianProduct(H,E) };
P := PerpSets(SympSp3);
perpsetsPointsIntersections := 
{ elt[1][2] meet elt[2][2] : elt in CartesianProduct(P,P) 
  | InnerProduct(elt[1][1],elt[2][1]) eq 1 };
lines := QuantumSubspaces(SympSp3,1);
perpsets3LinesIntersections := 
{ 
  { line : line in lines | line subset geometryPoints } : 
  geometryPoints in perpsetsPointsIntersections 
};
Doilies3 := perpsets3LinesIntersections join quadrics3Intersections;
"Number of doilies in W(5,2):";
#Doilies3;
increased3Doilies := 
&join
{ 
  &join
  {
    &join
    { 
      {
        IncreaseGeometry(geometry,{<pos,op,hyperplane>}) 
          : pos in [1..n]
      } : op in PauliOperators
    } : hyperplane in InnerOvoids(geometry) join InnerGrids(geometry) join InnerPerpsets(geometry)
  } : geometry in Doilies3
};
"Number of doilies of W(7,2) resulting from the doilies of W(5,2):";
#increased3Doilies;
"Number of doilies of W(7,2) resulting from the doilies of W(5,2) and also
recoverable by the doily of W(3,2):";
#(increased3Doilies meet increased2Doilies);
Doilies4 := quadricsIntersections join perpsetsLinesIntersections; // these are not doilies !
// todo, check that they are doilies
"Number of genuine W(7,2) doilies:";
#((Doilies4 diff increased2Doilies) diff increased3Doilies);
"Total number of doilies in W(7,2):";
#Doilies4; // this number can be calculated in doubt

// -------------------------------------------------------------------------- //
"
 --- Study of points order per signature --- ";
allSignatures :=
{ 
  <geometry,SHSignature(geometry)> 
    : geometry in perpsetsLinesIntersections join quadricsIntersections
};
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

exit;