freeze;

intrinsic QuantumSymplecticSpace(n::RngIntElt) -> ModTupRng
{ Builds the Symplectic space W(2*n-1,2) where n is the number of qubits 
  considered. }
  requirege n, 1;
  F2 := FiniteField(2);
  J := Matrix(F2, 2, 2, 
   [0,1,
    1,0]);
  In := ScalarMatrix(n, F2.1);
  Jn := TensorProduct(In,J);
  return SymplecticSpace(Jn);
end intrinsic;

intrinsic QSS(n::RngIntElt) -> ModTupRng
{ Short for QuantumSymplecticSpace. }
  return QuantumSymplecticSpace(n);
end intrinsic;

intrinsic QuantumInc(SympSp::ModTupFld) -> Inc
{ Computes an incidence structure from the symplectic space *SympSp*, whose elements are 
  the index of the elements of *SympSp* and whose blocks are the cliques of
  its commutation graph, and thus the generators of the symplectic space. }
  elts := {@ elt : elt in SympSp @};
  eltsIndex := {@ i : i in [1..#elts] | not IsZero(elts[i]) @};
  edgesIndex := { { ij[1], ij[2] } : ij in car<eltsIndex,eltsIndex> | 
    IsZero(InnerProduct(elts[ij[1]],elts[ij[2]])) and ij[1] lt ij[2] };
  graph := Graph<eltsIndex | edgesIndex>;
  cliquesIndex := [{eltsIndex[Index(vertex)] : vertex in clique} : 
                    clique in AllCliques(graph)];
  return IncidenceStructure<eltsIndex | cliquesIndex>;
end intrinsic;

intrinsic QuantumInc(n::RngIntElt) -> Inc
{ Computes incidence structures from a symplectic space, where the points are 
  the index of the elements for SympSp. The symplectic space itself is fixed
  by the number of qubits n. }
  return QuantumInc(QuantumSymplecticSpace(n));
end intrinsic;

intrinsic PauliOperatorInterpretation(vect::ModTupFldElt,interpretation::SeqEnum[FldCycElt[FldRat]]) 
  -> AlgMat
{ Returns the Pauli operator corresponding to the given vector form a 
  symplectic space for the interpretation given. The phases correspond to the 
  base matrices in the corresponding order: I,Z,X,Y. The canonical 
  interpretation is [K!1,1,1,-i]. }
  require #interpretation eq 4: "The interpretation must be the phase \
associated to each combination of 0s and 1s, with the resulting operator being \
on each qubit interpretation_j*Z^{a_i}*X^{b_i} and j corresponding to 1 for \
(0,0), 2 for (0,1), 3 for (1,0) and 4 for (1,1)";
  n := Ncols(vect) div 2;
  Int := IntegerRing();
  K<i> := Parent(interpretation[1]); 
  require &and{ element^4 eq 1 : element in interpretation}: "The elements of \
the interpretation must be in {1,-1,i,-i}";
  result := ScalarMatrix(K, 1, 1);
  X := Matrix(K, 2, 2, [0,1,1,0]);
  Z := DiagonalMatrix(K, 2, [1,-1]);
  for index := 1 to n do
    a := Int!vect[2*index-1];
    b := Int!vect[2*index];
    j := a*1+b*2+1;
    result := TensorProduct(result,interpretation[j]*Z^a*X^b);
  end for;
  return result;
end intrinsic;

intrinsic PauliOperatorPower(vect::ModTupFldElt) -> AlgMat
{ Returns the Pauli operator corresponding to the given vector from a 
  symplectic space.
  PauliOperatorPower(a)=PauliOperatorInterpresation(a,[1,1,1,1]) }
  n := Ncols(vect) div 2;
  Int := IntegerRing();
  K<i> := CyclotomicField(4);
  result := ScalarMatrix(K, 1, 1);
  X := Matrix(K, 2, 2, [0,1,1,0]);
  Z := DiagonalMatrix(K, 2, [1,-1]);
  for index := 1 to n do
    result := TensorProduct(result,Z^Int!vect[2*index-1]*X^Int!vect[2*index]);
  end for;
  return result;
end intrinsic;

intrinsic PauliOperatorCanonical(vect::ModTupFldElt) -> AlgMat
{ Returns the Pauli operator corresponding to the given vector from a 
  symplectic space. 
  PauliOperatorCanonical(a)=PauliOperatorInterpresation(a,[1,1,1,-i]) }
  n := Ncols(vect) div 2;
  Int := IntegerRing();
  K<i> := CyclotomicField(4);
  result := ScalarMatrix(K, 1, 1);
  I := DiagonalMatrix(K, 2, [1,1]);
  X := Matrix(K, 2, 2, [0,1,1,0]);
  Y := Matrix(K, 2, 2, [0,-i,i,0]);
  Z := DiagonalMatrix(K, 2, [1,-1]);
  // operator := I;
  for index := 1 to n do
    if [vect[2*index-1], vect[2*index]] eq [1,1] then
      operator := Y;
    else
      operator := Z^Int!vect[2*index-1]*X^Int!vect[2*index];
    end if;
    result := TensorProduct(result,operator);
  end for;
  return result;
end intrinsic;

intrinsic AntiCommute(v::ModTupFldElt,w::ModTupFldElt) -> BoolElt
{ }
  A := PauliOperatorCanonical(v);
  B := PauliOperatorCanonical(w);
  return A*B eq -B*A;
end intrinsic;

intrinsic SetSign(s::SetEnum[ModTupFldElt],f::UserProgram) -> RngIntElt
{ Computes the sign of a set of elements of a quantum symplectic space given the
  interpretation function f. }
  require &+s eq Parent(Random(s))!0: "The set must be linearly dependent.";
  matProd := &*{f(elt) : elt in s};
  require #Eigenvalues(matProd) eq 1: 
    "Error in computation, there must be only one eigenvalue in set product, " 
      cat Sprint(#Eigenvalues(matProd)) cat " found";
  return Random(Eigenvalues(matProd))[1];
end intrinsic;

intrinsic SetSign(s::SetEnum[ModTupFldElt],f::Intrinsic) -> RngIntElt
{ Computes the sign of a set of elements of a quantum symplectic space given the
  interpretation intrinsic f. }
  function fRelaxed(p)
    return f(p);
  end function;
  return SetSign(s,fRelaxed);
end intrinsic;

intrinsic SetSign(s::SetEnum[ModTupFldElt]) -> RngIntElt
{ Computes the sign of a set of elements of a quantum symplectic space for the
  canonical interpretation. }
  return SetSign(s,PauliOperatorCanonical);
end intrinsic;

intrinsic SetSign(s::SetIndx[ModTupFldElt]) -> RngIntElt
{ Computes the sign of a set of elements of a quantum symplectic space for the
  canonical interpretation. }
  return SetSign(IndexedSetToSet(s),PauliOperatorCanonical);
end intrinsic;

// the usage of the blocks is meant to avoid commutation checks between every
// pairs of points, but for low dimensions, it is less efficient because of
// the big number of blocks and the low number of points in each subspace, this
// needs to be further investigated.
intrinsic QuantumSubspaces(SympSp::ModTupFld,k::RngIntElt) -> SetEnum[SetEnum[ModTupFldElt]]
{ Computes the set of all isotropic subspaces of dimension k for the 
  symplectic space SympSp. }
  SympInc := QuantumInc(SympSp);
  points := {@ elt : elt in SympSp @};
  posPoints := { elt : elt in SympSp | not IsZero(elt) };
  B := Blocks(SympInc);
  n := Degree(SympSp)/2;
  // recursion initialization and special cases optimizations   
  if k eq 0 then 
    return { {point} : point in points | not IsZero(point) };
  elif k eq 1 then
    return { pair join {&+pair} : pair in Subsets(posPoints,2) |
      InnerProduct(Setseq(pair)[1],Setseq(pair)[2]) eq 0 };
      // there may be a problem here, the result of Setseq is not guarantied to
      // be the same at each call ?
  elif k eq n-1 then
    return { { points[i] : i in Support(block) } : block in B };
  end if;
  subspaces := {};
  previousSubspace := QuantumSubspaces(SympSp,k-1);
  for block in B do
    for blockSubsetIndex in Subsets(Support(block),k+1) do
      blockSubset := {points[index] : index in blockSubsetIndex};
      if &and{ not blockSubset subset elt : elt in previousSubspace } then 
        subspaces join:= { Closure(blockSubset,'+') diff { SympSp!0 } };
      end if;
    end for;
  end for;
  return subspaces;
end intrinsic;

intrinsic Closure(s::SetEnum,f::Intrinsic) -> SetEnum
{ Computes the closure of s by f. }
  previousStep := {};
  result := s;
  while previousStep ne result do
    previousStep := result;
    for pair in Subsequences(previousStep,2) do
      result join:= { f(pair[1],pair[2]) };
    end for;
  end while;
  return result;
end intrinsic;

intrinsic Lines(SympSp::ModTupFld) -> SetEnum[SetEnum[ModTupFldElt]]
{ Computes the set of all isotropic subspaces of dimension 1 for the 
  symplectic space SympSp. }
  return QuantumSubspaces(SympSp,1);
end intrinsic;

intrinsic Generators(SympSp::ModTupFld) -> SetEnum[SetEnum[ModTupFldElt]]
{ Computes the set of all isotropic subspaces of dimension n-1 for the 
  symplectic space SympSp. }
  n := Degree(SympSp)/2;
  return QuantumSubspaces(SympSp,n-1);
end intrinsic;

intrinsic QuantumSubspaces(n::RngIntElt,k::RngIntElt) -> SetEnum[SetEnum[ModTupFldElt]]
{ Computes the list of isotropic subspaces of dimension k for a system of n  
  qubits. }
  requirerange k, 0, n-1;
  SympSp := QuantumSymplecticSpace(n);
  return QuantumSubspaces(SympSp,k);
end intrinsic;

intrinsic SymplecticPoint(n::RngIntElt,p::RngIntElt) -> ModTupFldElt
{ Returns the vector of the symplectic space corresponding to the pth point of 
  the quantum symplectic space on n qubits. }
  return [ elt : elt in QuantumSymplecticSpace(n) ][p];
end intrinsic;

intrinsic SymplecticPoint(SympSp::ModTupFld,p::RngIntElt) -> ModTupFldElt
{ Returns the vector of the symplectic space corresponding to the pth point of 
  the quantum symplectic space SympSp. }
  return [ elt : elt in SympSp ][p];
end intrinsic;

intrinsic SymplecticIndex(p::ModTupFldElt) -> RngIntElt
{ Returns the index of the vector p in the list of its containing symplectic 
  space (the vector must be from a finite space). }
  return Index([ elt : elt in Parent(p) ], p);
end intrinsic;

intrinsic QuantumSubspacesIndex(n::RngIntElt,k::RngIntElt) -> SetEnum[SetEnum[RngIntElt]]
{ Computes the list of isotropic subspaces of dimension k for a system of n  
  qubits. Instead of returning the points of the symplectic space, return their 
  index for a more compact representation. }
  requirerange k, 0, n-1;
  SympSp := QuantumSymplecticSpace(n);
  SympSpElts := [ elt : elt in SympSp ];
  Subspaces := QuantumSubspaces(SympSp,k);
  return {{Index(SympSpElts, elt) : elt in subspace} : subspace in Subspaces};
end intrinsic;

intrinsic ToPauliString(p::ModTupFldElt) -> MonStgElt
{ Translates the element of the symplectic space to a string composed of the 
  corresponding Pauli matrices }
  pSeq := Eltseq(p);
  str := "";
  for index in [1..#pSeq/2] do
    if pSeq[2*index-1..2*index] eq [0,0] then
      str *:= "I";
    elif pSeq[2*index-1..2*index] eq [0,1] then
      str *:= "X";
    elif pSeq[2*index-1..2*index] eq [1,0] then
      str *:= "Z";
    elif pSeq[2*index-1..2*index] eq [1,1] then
      str *:= "Y";
    end if;
  end for;
  return str;
end intrinsic;

intrinsic TPS(p::ModTupFldElt) -> MonStgElt
{ Short for ToPauliString. }
  return ToPauliString(p);
end intrinsic;

intrinsic FromPauliString(str::MonStgElt) -> ModTupFldElt
{ Creates a point of the symplectic space from a string describing an operator 
  of the Pauli group. }
  require &and{ char in "XYZI" : char in Eltseq(str) }: "String must be \
composed of the Pauli operators (X, Y, Z and I)" ;
  p := [];
  for index in [1..#str] do
    if str[index] eq "I" then
      p cat:= [0,0];
    elif str[index] eq "X" then
      p cat:= [0,1];
    elif str[index] eq "Y" then
      p cat:= [1,1];
    elif str[index] eq "Z" then
      p cat:= [1,0];
    end if;
  end for;
  return QSS(#str)!p;
end intrinsic;

intrinsic FPS(str::MonStgElt) -> ModTupFldElt
{ Short for FromPauliString. }
  return FromPauliString(str);
end intrinsic;

intrinsic PrettyPrintSigned(geometry::SetEnum[SetEnum[ModTupFldElt]])
  -> MonStgElt
{ Display the geometry using the Pauli matrices representatives I, X, Y and Z 
  aswell as the sign of each line given the standard interpretation. }
  str := "{\n";
  iteration := 1;
  for l in geometry do
    str := Sprintf("    %o%o -> %2o\n",{ ToPauliString(p) : p in l },
      iteration eq #geometry select " " else ",", SetSign(l)); 
    iteration +:= 1;
  end for;
  str *:= "}";
  return str;
end intrinsic;

intrinsic PrettyPrint(geometry::SetEnum[SetEnum[ModTupFldElt]])
  -> SetEnum[SetEnum[MonStgElt]]
{ Display the geometry using the Pauli matrices representatives I, X, Y and Z. }
  return { PrettyPrint(l): l in geometry };
end intrinsic;

intrinsic PrettyPrint(context::SetEnum[ModTupFldElt]) -> SetEnum[MonStgElt]
{ Display the geometry using the Pauli matrices representatives I, X, Y and Z. }
  return { ToPauliString(p) : p in context };
end intrinsic;