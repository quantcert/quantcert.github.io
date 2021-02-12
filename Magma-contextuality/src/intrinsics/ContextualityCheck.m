freeze;

intrinsic IsContextual(geometry::SetEnum[SetEnum[ModTupFldElt]]) -> BoolElt
{ Checks contextuality using the reformulation of the problem in terms of linear
  algebra}
  Z2 := GF(2);
  points := SetToIndexedSet(&join geometry);
  inc := IncidenceStructure<points|geometry>;
  mat := Matrix(Z2, IncidenceMatrix(inc));
  vect := Vector(Z2, [ (SetSign(subspace,PauliOperatorCanonical) eq -1) 
    select 1 else 0 : subspace in geometry ]);
  // Note : we are shifting from ({-1,1},x,1) to ({Z/2Z},+,0) for the "sign" of
  // the lines, hence the fact that -1 maps to 0 
  isConsistent, solution, kernel := IsConsistent(mat,vect);
  return not isConsistent;
end intrinsic;
