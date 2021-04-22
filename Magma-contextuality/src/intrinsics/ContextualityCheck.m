freeze;

intrinsic IsContextual(geometry::SetEnum[SetEnum[ModTupFldElt]],interpretation::UserProgram) 
  -> BoolElt
{ Checks contextuality using the reformulation of the problem in terms of linear
  algebra, given the interpretation function (the interpretation function must 
  transform point of the quantum symplectic space into a element of the Pauli 
  group. }
  Z2 := GF(2);
  points := SetToIndexedSet(&join geometry);
  inc := IncidenceStructure<points|geometry>;
  mat := Matrix(Z2, IncidenceMatrix(inc));
  vect := Vector(Z2, [ (SetSign(subspace,interpretation) eq -1) 
    select 1 else 0 : subspace in geometry ]);
  isConsistent, solution, kernel := IsConsistent(mat,vect);
  return not isConsistent;
end intrinsic;

intrinsic IsContextual(geometry::SetEnum[SetEnum[ModTupFldElt]],interpretation::Intrinsic) 
  -> BoolElt
{ Checks contextuality using the reformulation of the problem in terms of linear
  algebra, given the interpretation function (the interpretation function must 
  transform point of the quantum symplectic space into a element of the Pauli 
  group. }
  function interpretationRelaxed(p)
    return interpretation(p);
  end function;
  return IsContextual(geometry,interpretationRelaxed);
end intrinsic;

intrinsic IsContextual(geometry::SetEnum[SetEnum[ModTupFldElt]],interpretation::SeqEnum[FldCycElt[FldRat]]) 
  -> BoolElt
{ Checks contextuality using the reformulation of the problem in terms of linear
  algebra, given the interpretation. The interpretation in this case must be 
  a list of four elements, in [1,-1,i,-i] corresponding to the second argument 
  of PauliOperatorInterpretation. }
  function PauliOperator(p)
    return PauliOperatorInterpretation(p,interpretation);
  end function;
  return  IsContextual(geometry,PauliOperator);
end intrinsic;

intrinsic IsContextual(geometry::SetEnum[SetEnum[ModTupFldElt]]) -> BoolElt
{ Checks contextuality using the reformulation of the problem in terms of linear
  algebra. }
  return IsContextual(geometry,PauliOperatorCanonical);
end intrinsic;

intrinsic NumberOfI(point::ModTupFldElt) -> RngIntElt
{ Return the number of Identities in a point, e.g. 1 for XIX and 2 for IIZ .}
  SympSp := Parent((point));
  n := Degree(SympSp)/2;
  return #{ i : i in [1..n] | point[2*i-1] eq 0 and point[2*i] eq 0 };
end intrinsic;

intrinsic SHSignature(geometry::SetEnum[SetEnum[ModTupFldElt]]) -> Tup
{ The signature *(C_i,Oa,Ob,Oc)* of *geometry* as defined by Metod and Holweck. 
  - C_i is the number of negative lines of the doily;
  - Oa is the number of observables of the doily which contain two identity 
    (example IIX or IYI or XII....);
  - Ob is the number of observables of the doily which contain one identity 
    (example XIX, ZYI,...);
  - Oc is the number of observables of the doily which contain no identity 
    (XXX, XYZ,....). }
  points := &join geometry;
  SympSp := Parent(Random(points));
  n := Degree(SympSp)/2;
  C_i := #{ line : line in geometry | SetSign(line) eq -1 };
  function Oi(pointSet, nbId)
    return #{ point : point in pointSet | NumberOfI(point) eq nbId };
  end function;
  return <C_i,[ Oi(points,k) : k in [ n-1 .. 0 by -1 ] ]>;
end intrinsic;

intrinsic SubspaceComparison(ssp1::SetIndx[ModTupFldElt],ssp2::SetIndx[ModTupFldElt]) -> RngIntElt
{ Compare two subspaces so they can be consistently sorted. }
  index := 1;
  while index lt Degree(Random(ssp1)) do
    if ssp1[index] lt ssp2[index] then 
      return -1;
    elif ssp1[index] gt ssp2[index] then
      return 1;
    end if;
    index +:= 1;
  end while;
  return 0;
end intrinsic;

intrinsic SubspaceComparison(ssp1::SetEnum[ModTupFldElt],ssp2::SetEnum[ModTupFldElt]) -> RngIntElt
{ Compare two subspaces so they can be consistently sorted. }
  return SubspaceComparison(SetToIndexedSet(ssp1),SetToIndexedSet(ssp2));
end intrinsic;

intrinsic ContextualityDegree(geometry::SetEnum[SetEnum[ModTupFldElt]] : debug:=false) -> RngIntElt
{ Returns the contextuality degree of a geometry where the contextuality degree 
  is the minimum number of constraints non solvable for the system *Ax=v*. In 
  this system, *A* is the incidence matrix of the geometry, *v* is the value 
  associated to each context, and *x* is the unknown. }

  subspaces := 
    Sort(
      {@ Sort(SetToIndexedSet(subsp)) : subsp in geometry @},
      func<a,b|SubspaceComparison(a,b)>
    );
  points := Sort(SetToIndexedSet(&join subspaces));
  inc := IncidenceStructure<points|subspaces>;
  Z2 := GF(2);
  mat := Matrix(Z2, IncidenceMatrix(inc));
  vect := Vector(
    Z2, 
    [ (SetSign(subspace) eq -1) select 1 else 0 : subspace in subspaces ]
  );

  return  Min({ Distance(vect,Parent(vect)!elt) : elt in Image(mat) });
end intrinsic;

intrinsic ContextualityDegree(geometry::SetIndx[SetIndx[ModTupFldElt]] : debug:=false) -> RngIntElt
{ Returns the contextuality degree of a geometry where the contextuality degree 
  is the minimum number of constraints non solvable for the system *Ax=v*. In 
  this system, *A* is the incidence matrix of the geometry, *v* is the value 
  associated to each context, and *x* is the unknown. }
  return ContextualityDegree({ IndexedSetToSet(subsp) : subsp in geometry });
end intrinsic;
