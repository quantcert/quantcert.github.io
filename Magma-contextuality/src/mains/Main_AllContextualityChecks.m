geometries := [<WLines,"WL">, <WBlocks,"WB">, <PerpSets,"PS">, <Elliptics,"El">, 
  <Hyperbolics,"Hy">];

"
Data is provided in the following format: 
  |family(size):contextuality(number of instances)|time|
For a given size, INC means the contextuality is inconsistent in the family, 
CON means all members of the family are contextual and NC means no member of the
family is contextual. The time given is the time to compute the results for the 
whole family in seconds.
Families are shorthanded as follows:
WLines:WL, WBlocks:WB, PerpSets:PS, Elliptics:El, Hyperbolics:Hy

|contextuality of famillies|time(s)|
|--------------------------|-------|";

for size := 2 to 5 do
  SympSp := QuantumSymplecticSpace(size);

  for geometry in geometries do
    if not (geometry[1] eq Elliptics and size eq 2) then
      t := Realtime();
      contextuality := [ IsContextual(instance[3]) : 
        instance in geometry[1](SympSp) ];
      
      allIdentical := #SequenceToSet(contextuality) eq 1;
      
      Sprintf("|%o(%o):%o(x%o)|%o|",geometry[2], size,
        allIdentical select contextuality[1] select "CON" else "NC " else "INC", 
        #contextuality, Realtime(t)
      );
    end if;
  end for;

end for;

exit;
