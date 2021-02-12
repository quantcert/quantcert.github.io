intrinsic Quadric(SympSp::ModTupFld, QuadForm::UserProgram, lines::SetEnum[SetEnum[ModTupFldElt]]) 
  -> SetEnum[ModTupFldElt], SetEnum[SetEnum[ModTupFldElt]]
{ Computes the plane defined by the quadratic form *QuadForm*. This plane is 
  also called a quadric. The set of those points is returned as well as
  the geometry formed by the lines composed of those points. }
  quadricPoints := {a : a in SympSp | IsZero(QuadForm(a)) and not IsZero(a)};
  return quadricPoints, { line : line in lines | line subset quadricPoints };
end intrinsic;

intrinsic Quadric(point) -> SetEnum[ModTupFldElt], SetEnum[SetEnum[ModTupFldElt]]
{ Computes the quadric defined by the point *point*. The set of those points is 
  returned as well as the geometry formed by the lines composed of those 
  points. }
  SympSp := Parent(point);
  QuadForm := QuadraticForm(point);
  lines := IsotropicSubspaces(SympSp,1);
  return Quadric(SympSp, QuadForm, lines);
end intrinsic;

intrinsic BaseQuadraticForm(SympSp::ModTupFld) -> UserProgram
{ Computes the quadratic form }
  n := Degree(SympSp);
  require IsEven(n): "*SympSp* must be a QuatumSymplecticSpace, and hence have \
a pair degree";
  function Q(x)
    return &+[x[2*i-1]*x[2*i] : i in {1 .. n/2}];
  end function;
  return  Q;
end intrinsic;

intrinsic QuadraticForm(point::ModTupFldElt) -> UserProgram
{ Computes the quadratic form null at the point *point*. }
  SympSp := Parent(point);
  Q := BaseQuadraticForm(SympSp);
  function Qp(x)
    return  Q(x) + InnerProduct(x,point);
  end function;
  return  Qp;
end intrinsic;

intrinsic Hyperbolics(SympSp::ModTupFld) -> SeqEnum
{ Computes the hyperbolics in SympSp. Returns a list of tuples, where each 
  tuples is a point generating the hyperboloid, the points of this hyperboloid
  and the lines of this hyperboloid. }
  BaseQuadricPoints, BaseQuadric := Quadric(SympSp!0);
  GQ := [<SympSp!0, BaseQuadricPoints, BaseQuadric>];
  for point in BaseQuadricPoints do
    qPoints, qLines := Quadric(point);
    Append(~GQ,<point, qPoints, qLines>);
  end for;
  return  GQ;
end intrinsic;

intrinsic Elliptics(SympSp::ModTupFld) -> SeqEnum
{ Computes the elliptics in SympSp. Returns a list of tuples, where each 
  tuples is a point generating the ellipsoid, the points of this ellipsoid
  and the lines of this ellipsoid. }
  require Degree(SympSp)/2 gt 2 : "There is no Elliptic for less than 3 qubits\
 (the given space has " cat Sprint(Degree(SympSp)/2) cat " qubits)"; 
  GQ := [];
  BaseQuadricPoints, BaseQuadric := Quadric(SympSp!0);
  for point in {point : point in SympSp } diff (BaseQuadricPoints join {SympSp!0}) do
    qPoints, qLines := Quadric(point);
    Append(~GQ,<point, qPoints, qLines>);
  end for;
  return  GQ;
end intrinsic;

intrinsic PerpSet(point::ModTupFldElt) -> SetEnum[SetEnum[ModTupFldElt]]
{ Return the set of all lines containing the point *point*. This is called 
  PerpSet, or the plan concurrent to *point*. }
  require not IsZero(point):"point cannot be null";
  SympSp := Parent(point);
  lines := IsotropicSubspaces(SympSp,1);
  return { line : line in lines | point in line };
end intrinsic;

intrinsic PerpSets(SympSp::ModTupFld) -> SeqEnum
{ Computes the PerpSets in SympSp. Returns a list of tuples, where each 
  tuples is a point generating the PerpSet, the points of this PerpSet
  and the lines of this PerpSet. }
  P := [];
  for point in { elt : elt in SympSp | not IsZero(elt) } do
    ps := PerpSet(point);
    psPoints := &join ps;
    Append(~P,<point, psPoints, ps>);
  end for;
  return P;
end intrinsic;

intrinsic WLines(SympSp::ModTupFld) -> SeqEnum
{ Computes W(2*n-1,2) with the geometry being its lines. Returns a list of 
  tuples (event though this list contains a single element, to be coherent with 
  the other methods of this file), the point has no particular signification, it 
  is arbitrarily chosen to be the null vector.. }
  points := { elt : elt in SympSp | not IsZero(elt) };
  lines := IsotropicSubspaces(SympSp,1);
  return [<SympSp!0,points,lines>];
end intrinsic;

intrinsic WBlocks(SympSp::ModTupFld) -> SeqEnum
{ Computes W(2*n-1,2) with the geometry being its lines. Returns a list of 
  tuples (event though this list contains a single element, to be coherent with 
  the other methods of this file), the point has no particular signification, it 
  is arbitrarily chosen to be the null vector.. }
  points := { elt : elt in SympSp | not IsZero(elt) };
  incBlocks := Blocks(QuantumInc(SympSp));
  blocks := { { SymplecticPoint(SympSp,i) : i in Support(block) } : block in incBlocks };
  return [<SympSp!0,points,blocks>];
end intrinsic;
