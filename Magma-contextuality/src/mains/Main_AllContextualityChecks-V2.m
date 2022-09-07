geometries := [<WLines,"WL">, <WBlocks,"WB">, <PerpSets,"PS">, <Elliptics,"El">, 
  <Hyperbolics,"Hy">];

"
Data is provided in the given format : 
  family(size):contextuality(number of instances)
For a given size, INC means the contextuality is inconsistent in the family, 
CON means all member of the family are contextual and NC means no member of the
family is contextual.
Families are shorthanded as follows:
WLines:WL, WBlocks:WB, PerpSets:PS, Elliptics:El, Hyperbolics:Hy
";

for size := 2 to 5 do
  SympSp := QuantumSymplecticSpace(size);

  for geometry in geometries do
    if not (geometry[1] eq Elliptics and size eq 2) then
      t := Realtime();
      contextuality := [ IsContextual(instance[3]) : 
        instance in geometry[1](SympSp) ];
      
      allIdentical := #SequenceToSet(contextuality) eq 1;
      
      if allIdentical then
        Sprintf("%o(%o):%o(x%o)",geometry[2], size,
          contextuality[1] select "CON" else "NC ", #contextuality);
      else
        Sprintf("%o(%o):INC(x%o)",geometry[2], size, #contextuality);
      end if;
      Realtime(t);
    end if;
  end for;

end for;

exit;
