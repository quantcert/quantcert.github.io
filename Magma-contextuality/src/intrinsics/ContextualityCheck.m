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

intrinsic IsContextual(geometry::SetEnum[SetEnum[ModTupFldElt]]) -> BoolElt
{ Checks contextuality using the reformulation of the problem in terms of linear
  algebra. }
  return IsContextual(geometry,PauliOperatorCanonical);
end intrinsic;