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

intrinsic QuantumInc(SympSp::ModTupFld) -> Inc
{ Computes incidence structures from a symplectic space, where the points are 
  the index of the elements for SympSp. }
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
  for i := 1 to n do
    if [vect[2*i-1], vect[2*i]] eq [1,1] then
      operator := Y;
    else
      operator := Z^Int!vect[2*i-1]*X^Int!vect[2*i];
    end if;
    result := TensorProduct(result,operator);
  end for;
  return result;
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
  require &+s eq Parent(Random(s))!0: "The set must be linearly dependent.";
  matProd := &*{f(elt) : elt in s};
  require #Eigenvalues(matProd) eq 1: 
    "Error in computation, there must be only one eigenvalue in set product, " 
      cat Sprint(#Eigenvalues(matProd)) cat " found";
  return Random(Eigenvalues(matProd))[1];
end intrinsic;

intrinsic IsotropicSubspaces(SympSp::ModTupFld,k::RngIntElt) -> SetEnum[SetEnum[ModTupFldElt]]
{ Computes the set of all isotropic subspaces of dimension k for the 
  symplectic space SympSp. }
  SympInc := QuantumInc(SympSp);
  points := {@ elt : elt in SympSp @};
  B := Blocks(SympInc);
  if k eq 0 then 
    return { {point} : point in points | not IsZero(point) };
  end if;
  subspaces := {};
  previousSubspace := IsotropicSubspaces(SympSp,k-1);
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

intrinsic IsotropicSubspaces(n::RngIntElt,k::RngIntElt) -> SetEnum[SetEnum[ModTupFldElt]]
{ Computes the list of isotropic subspaces of dimension k for a system of n  
  qubits. }
  requirerange k, 0, n-1;
  SympSp := QuantumSymplecticSpace(n);
  return IsotropicSubspaces(SympSp,k);
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

intrinsic IsotropicSubspacesIndex(n::RngIntElt,k::RngIntElt) -> SetEnum[SetEnum[RngIntElt]]
{ Computes the list of isotropic subspaces of dimension k for a system of n  
  qubits. Instead of returning the points of the symplectic space, return their 
  index for a more compact representation. }
  requirerange k, 0, n-1;
  SympSp := QuantumSymplecticSpace(n);
  SympSpElts := [ elt : elt in SympSp ];
  Subspaces := IsotropicSubspaces(SympSp,k);
  return {{Index(SympSpElts, elt) : elt in subspace} : subspace in Subspaces};
end intrinsic;
