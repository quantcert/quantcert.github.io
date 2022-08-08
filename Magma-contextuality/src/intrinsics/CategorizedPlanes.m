freeze;

intrinsic Quadric(SympSp::ModTupFld, QuadForm::UserProgram, lines::SetEnum[SetEnum[ModTupFldElt]]) 
  -> SetEnum[ModTupFldElt], SetEnum[SetEnum[ModTupFldElt]]
{ Computes the plane defined by the quadratic form *QuadForm*. This plane is 
  also called a quadric. The set of those points is returned as well as
  the geometry formed by the lines composed of those points. }
  quadricPoints := {a : a in SympSp | IsZero(QuadForm(a)) and not IsZero(a)};
  return quadricPoints, { line : line in lines | line subset quadricPoints };
end intrinsic;

intrinsic Quadric(point::ModTupFldElt, lines::SetEnum[SetEnum[ModTupFldElt]]) 
  -> SetEnum[ModTupFldElt], SetEnum[SetEnum[ModTupFldElt]]
{ Computes the quadric defined by the point *point*. The set of those points is 
  returned as well as the geometry formed by the lines composed of those 
  points. }
  SympSp := Parent(point);
  QuadForm := QuadraticForm(point);
  return Quadric(SympSp, QuadForm, lines);
end intrinsic;

intrinsic Quadric(point::ModTupFldElt) 
  -> SetEnum[ModTupFldElt], SetEnum[SetEnum[ModTupFldElt]]
{ Computes the quadric defined by the point *point*. The set of those points is 
  returned as well as the geometry formed by the lines composed of those 
  points. Short version of *Quadric(point,lines)*, if called multiple times, the 
  lines will be recomputed each time, so avoid it if possible.}
  lines := IsotropicSubspaces(Parent(point),1);
  return Quadric(point, lines);
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

intrinsic Hyperbolics(SympSp::ModTupFld) -> SeqEnum[Tup]
{ Computes the hyperbolics in SympSp. Returns a list of tuples, where each 
  tuples is a point generating the hyperbolic, the points of this hyperbolic
  and the lines of this hyperbolic. }
  BaseQuadricPoints, BaseQuadric := Quadric(SympSp!0);
  GQ := {<SympSp!0, BaseQuadricPoints, BaseQuadric>};
  lines := IsotropicSubspaces(SympSp,1);
  for point in BaseQuadricPoints do
    qPoints, qLines := Quadric(point,lines);
    GQ join:= {<point, qPoints, qLines>};
  end for;
  return GQ;
end intrinsic;

intrinsic Elliptics(SympSp::ModTupFld) -> SetEnum[Tup]
{ Computes the elliptics in SympSp. Returns a list of tuples, where each 
  tuples is a point generating the elliptic, the points of this elliptic
  and the lines of this elliptic. 
  Caution ! On two qubits, there is no line in an elliptic. }
  GQ := {};
  BaseQuadricPoints, BaseQuadric := Quadric(SympSp!0);
  lines := IsotropicSubspaces(SympSp,1);
  for point in { point : point in SympSp } diff (BaseQuadricPoints join {SympSp!0}) do
    qPoints, qLines := Quadric(point,lines);
    GQ join:= {<point, qPoints, qLines>};
  end for;
  return GQ;
end intrinsic;

intrinsic Quadrics(SympSp::ModTupFld) -> SetEnum[Tup], SetEnum[Tup]
{ Computes all the quadrics in SympSp. Returns two list of tuples, where each 
  tuples is a point generating the quadric, the points of this quadric and the 
  lines of this quadric. The first list contains all the hyperbolics and the 
  second list contains all the }
  Q0 := BaseQuadraticForm(SympSp);  
  BaseQuadricPoints := {a : a in SympSp | IsZero(Q0(a))};
  lines := IsotropicSubspaces(SympSp,1);
  hyperbolics := {};
  elliptics := {};
  for point in SympSp do
    qPoints, qLines := Quadric(point,lines);
    if point in BaseQuadricPoints then
      hyperbolics join:= {<point, qPoints, qLines>};
    else
      elliptics join:= {<point, qPoints, qLines>};
    end if;
  end for;
  return hyperbolics, elliptics;
end intrinsic;

intrinsic PerpSet(point::ModTupFldElt, lines::SetEnum[SetEnum[ModTupFldElt]]) 
  -> SetEnum[SetEnum[ModTupFldElt]]
{ Return the set of all lines containing the point *point*. This is called 
  PerpSet, or the plan concurrent to *point*. }
  require not IsZero(point):"point cannot be null";
  return { line : line in lines | point in line };
end intrinsic;

intrinsic PerpSet(point::ModTupFldElt) -> SetEnum[SetEnum[ModTupFldElt]]
{ Return the set of all lines containing the point *point*. This is called 
  PerpSet, or the plan concurrent to *point*. Short version of 
  *PerpSet(point,lines)*, if called multiple times, the lines will be recomputed
  each time, so avoid it if possible. }
  SympSp := Parent(point);
  lines := IsotropicSubspaces(SympSp,1);
  return PerpSet(point, lines);
end intrinsic;

intrinsic PerpSets(SympSp::ModTupFld) -> SeqEnum[Tup]
{ Computes the PerpSets in SympSp. Returns a list of tuples, where each 
  tuples is a point generating the PerpSet, the points of this PerpSet
  and the lines of this PerpSet. }
  P := {};
  lines := IsotropicSubspaces(SympSp,1);
  for point in { elt : elt in SympSp | not IsZero(elt) } do
    ps := PerpSet(point, lines);
    psPoints := &join ps;
    P join:= {<point, psPoints, ps>};
  end for;
  return P;
end intrinsic;

intrinsic WLines(SympSp::ModTupFld) -> SeqEnum[Tup]
{ Computes W(2*n-1,2) with the geometry being its lines. Returns a list of 
  tuples (even though this list contains a single element, to be coherent with 
  the other methods of this file), the point has no particular signification, it 
  is arbitrarily chosen to be the null vector.. }
  points := { elt : elt in SympSp | not IsZero(elt) };
  lines := IsotropicSubspaces(SympSp,1);
  return [<SympSp!0,points,lines>];
end intrinsic;

intrinsic WBlocks(SympSp::ModTupFld) -> SeqEnum[Tup]
{ Computes W(2*n-1,2) with the geometry being its lines. Returns a list of 
  tuples (even though this list contains a single element, to be coherent with 
  the other methods of this file), the point has no particular signification, it 
  is arbitrarily chosen to be the null vector.. }
  points := { elt : elt in SympSp | not IsZero(elt) };
  n := Integers()!(Degree(SympSp)/2);
  blocks := IsotropicSubspaces(SympSp, n-1);
  return [<SympSp!0,points,blocks>];
end intrinsic;

intrinsic IsDoily1(geometry::SetEnum[SetEnum[ModTupFldElt]]) -> BoolElt
{ Checks if the geometry correspond to the first definition of a doily :
  a doily has 15 points, 15 lines, each line has 3 points, each point is in 3 
  lines and all lines must have the product of their points equal to (0..0). }
  return
    #geometry eq 15
      and
    # &join geometry eq 15
      and
    &and{ #line eq 3 : line in geometry }
      and
    &and
    { 
      #{ line : line in geometry | point in line } eq 3 
        : point in &join geometry 
    } 
      and
    &and{ IsZero(&+ line) : line in geometry }
  ;
end intrinsic;