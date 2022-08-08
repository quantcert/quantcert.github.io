freeze;

intrinsic IncreasedPoint(finalSympSp::ModTupFld,p::ModTupFldElt,posObsNeutrals::SetEnum[Tup]) 
  -> ModTupFldElt
{ The point *p* with observables given by the second element of each pair of 
  *posObsNeutrals* as a list of two integers, added (as a tensor product on the
  decomposition) as position given by the second element of the corresponding 
  pair. }
  seqPoint := Eltseq(p);
  require #{ posObsNeutral[1] : posObsNeutral in posObsNeutrals } eq #posObsNeutrals : 
    "Multiple observable added at the same position";
  posObsPairsSeq := Setseq(posObsNeutrals);
  Sort(~posObsPairsSeq,func<a,b|a[1]-b[1]>); 
  // so the position to insert observables is already defined
  for posObsNeutral in posObsPairsSeq do
    pos := posObsNeutral[1];
    observable := p in posObsNeutral[3] select [0,0] else posObsNeutral[2];
    Insert(~seqPoint,2*pos-1,observable);
  end for;
  return finalSympSp!seqPoint;
end intrinsic;

intrinsic IncreasedPoint(p::ModTupFldElt,posObsNeutrals::SetEnum[Tup]) 
  -> ModTupFldElt
{ The point *p* with observables given by the second element of each pair of 
 Sort(~posObsPairsSeq,func<a,b|a[1]-b[1]>); *posObsNeutrals* as a list of two integers, added (as a tensor product on the
  decomposition) as position given by the second element of the corresponding 
  pair. }
  SympSp := Parent(p);
  n := Integers()!(Degree(SympSp)/2);
  finalSympSp := QuantumSymplecticSpace(n + #posObsNeutrals);
  return IncreasedPoint(finalSympSp,p,posObsNeutrals);
end intrinsic;

intrinsic IncreaseGeometry(geometry::SetEnum[SetEnum[ModTupFldElt]],posObsNeutrals::SetEnum[Tup]) 
  -> SetEnum[SetEnum[ModTupFldElt]]
{ The geometry with each point increased by the *IncreasedPoint* function. 
  *posObsNeutrals* is a set of pairs where the first element of the pair is the 
  position where to add the observable in second position. }
  SympSp := Parent(Random(Random(geometry)));
  n := Integers()!(Degree(SympSp)/2);
  finalSympSp := QuantumSymplecticSpace(n + #posObsNeutrals);
  return  
  { 
    { 
      IncreasedPoint(finalSympSp,p,posObsNeutrals) : p in context
    } : context in geometry 
  };  
end intrinsic;

intrinsic InnerOvoids(geometry::SetEnum[SetEnum[ModTupFldElt]]) 
  -> SetEnum[SetEnum[ModTupFldElt]]
{ The Ovoid included in the geometry. Careful, this is a very special case, it 
  probably works only for the doilies in dim 3... TODO generalize ? }
  ovoids := {};
  geometryPoints := &join geometry;
  for point in geometryPoints do
    nonCommutingPoints := { p : p in geometryPoints | InnerProduct(p,point) eq 1 } join { point } ;
    for ovoid in Subsets(nonCommutingPoints,5) do
      noneCommute := 
        &and{ InnerProduct(Setseq(pair)[1],Setseq(pair)[2]) eq 1
          : pair in Subsets(ovoid,2) };
      if noneCommute then
        ovoids join:= { ovoid };
        if #ovoids eq 6 then
          return ovoids;
        end if; 
      end if; 
    end for;
  end for;
  return  ovoids;
end intrinsic;

intrinsic InnerGrids(geometry::SetEnum[SetEnum[ModTupFldElt]]) 
  -> SetEnum[SetEnum[ModTupFldElt]]
{ The grid included in the geometry. Careful, this is a very special case, it 
  probably works only for the doilies in dim 3... TODO generalize ? }
  grids := {};
  concurrentLines := 
  &join
  { 
    { 
      { line1,line2 } : line1 in geometry | #(line1 meet line2) eq 1
    } : line2 in geometry 
  };
  for linesSet in concurrentLines do
    linesSeq := Setseq(linesSet);
    line1 := linesSeq[1];
    line2 := linesSeq[2];
    q := Random(line1 meet line2);
    for pair in CartesianProduct(line1,line2) do
      if q ne pair[1] and q ne pair[2] then
        p1 := pair[1];
        p2 := pair[2];
        p1b := Random(line1 diff {q,p1});
        p2b := Random(line2 diff {q,p2});
        for p3 in 
          { 
            p3 : p3 in &join geometry | 
            exists(p){p : p in &join geometry | {p,p1,p3} in geometry }
            and  
            exists(p){p : p in &join geometry | {p,p2,p3} in geometry }  
          } 
          do
          p4 := Random({ p : p in &join geometry | {p,p3,p1} in geometry });
          p5 := Random({ p : p in &join geometry | {p,p3,p2} in geometry });
          if exists(p){ p : p in &join geometry | {p4,p,p2b} in geometry }
            and
            exists(p){ p : p in &join geometry | {p5,p,p1b} in geometry }
          then
            p6 := Random({ p : p in &join geometry | {p4,p,p2b} in geometry });
            p7 := Random({ p : p in &join geometry | {p5,p,p1b} in geometry });
            grids join:= {{p1,p1b,q,p2,p2b,p3,p4,p5,p6,p7}};
          end if;
        end for;
      end if;
    end for;
  end for;
  return grids;
end intrinsic; 

intrinsic InnerPerpsets(geometry::SetEnum[SetEnum[ModTupFldElt]]) 
  -> SetEnum[SetEnum[ModTupFldElt]]
{ The perpsets included in the geometry. }
  return 
  { 
    { p : p in &join geometry | InnerProduct(p,point) eq 0 }
      : point in &join geometry 
  };
end intrinsic;

intrinsic Replace(seq::[],pos::RngIntElt,replacement::[]) -> []
{ Replaces the sequence starting at *pos* in *seq* by *replacement*. }
    return Insert(seq, pos, pos + #replacement - 1, replacement);
end intrinsic;

intrinsic Replace(~seq::[],pos::RngIntElt,replacement::[])
{ Replaces in place the sequence starting at *pos* in *seq* by *replacement*. }
    Insert(~seq, pos, pos + #replacement - 1, replacement);
end intrinsic;

intrinsic RelabelPoint(SympSp::ModTupFld,p::ModTupFldElt,pos::RngIntElt,newObservable::SeqEnum[RngIntElt]) 
  -> ModTupFldElt
{ Change the observable at position *pos* to be *newObservable*. }
  return  SympSp!Replace(Eltseq(p), 2*pos-1, newObservable);
end intrinsic;

intrinsic RelabelPoint(p::ModTupFldElt,pos::RngIntElt,newObservable::SeqEnum[RngIntElt]) 
  -> ModTupFldElt
{ Change the observable at position *pos* to be *newObservable*. }
  SympSp := Parent(p);
  return  RelabelPoint(SympSp,p,pos,newObservable);
end intrinsic;

intrinsic DeepPoints(geometry::SetEnum[SetEnum[ModTupFldElt]]) -> SetEnum[ModTupFldElt]
{ The set of deep points of a geometry, i.e. the points such as all the negative 
  lines of the geometry go through them. }
  return
  { 
    point 
      : point in &join geometry 
      | &and{ SetSign(line) eq -1 : line in geometry | point in line } 
  };
end intrinsic;

intrinsic ZeroPoints(geometry::SetEnum[SetEnum[ModTupFldElt]]) -> SetEnum[ModTupFldElt]
{ The set of deep points of a geometry, i.e. the points such as all the negative 
  lines of the geometry go through them. }
  return (&join geometry) diff (&join { line : line in geometry | SetSign(line) eq -1 });
end intrinsic;

intrinsic Order(point::ModTupFldElt,geometry::SetEnum[SetEnum[ModTupFldElt]],linesValues::SetEnum[Tup]) 
  -> RngIntElt
{ The order as the *point* of the *geometry*, as in, the number of negative 
  lines the point is in according to *lineValues*. }
  return 
  #{ 
    line 
      : line in geometry 
      | point in line 
          and 
        { lineValue[2] : lineValue in linesValues | lineValue[1] eq line } eq {-1} 
  };
end intrinsic;

intrinsic Order(point::ModTupFldElt,geometry::SetEnum[SetEnum[ModTupFldElt]],interpretation::Intrinsic) 
  -> RngIntElt
{ The order as the *point* of the *geometry*, as in, the number of negative 
  lines the point is in according to the *interpretation*. }
  linesValues := { <line,SetSign(line,interpretation)> : line in geometry };
  return Order(point,geometry,linesValues);
end intrinsic;

intrinsic Order(point::ModTupFldElt,geometry::SetEnum[SetEnum[ModTupFldElt]]) -> RngIntElt
{ The order as the *point* of the *geometry*, as in, the number of negative 
  lines the point is in. }
  return Order(point,geometry,PauliOperatorCanonical);
end intrinsic;

intrinsic PointsOrders(geometry::SetEnum[SetEnum[ModTupFldElt]],linesValues::SetEnum[Tup]) -> SetIndx[Tup]
{ Return the cardinal of each set to points corresponding to each order, as a 
  set of tuples where each tuple is composed of the order and the number of 
  points of this order in the geometry }
  pointsOrder := { <p,Order(p,geometry,linesValues)> : p in &join geometry };
  ordersCardinals :=
  {@
    < order, #{ pointOrder[1] : pointOrder in pointsOrder | pointOrder[2] eq order } >
      : order in [0..#geometry]
  @};  
  return {@ orderCardinal : orderCardinal in ordersCardinals | not orderCardinal[2] eq 0 @};
end intrinsic;

intrinsic PointsOrders(geometry::SetEnum[SetEnum[ModTupFldElt]],interpretation::Intrinsic) -> SetIndx[Tup]
{ Return the cardinal of each set to points corresponding to each order, as a 
  set of tuples where each tuple is composed of the order and the number of 
  points of this order in the geometry }
  linesValues := { <line,SetSign(line,interpretation)> : line in geometry };
  return PointsOrders(geometry,linesValues);
end intrinsic;

intrinsic PointsOrders(geometry::SetEnum[SetEnum[ModTupFldElt]]) -> SetIndx[Tup]
{ Return the cardinal of each set to points corresponding to each order, as a 
  set of tuples where each tuple is composed of the order and the number of 
  points of this order in the geometry }
  return PointsOrders(geometry,PauliOperatorCanonical);
end intrinsic;
