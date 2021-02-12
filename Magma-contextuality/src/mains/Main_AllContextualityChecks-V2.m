geometries := [<WLines,"WL">, <WBlocks,"WB">, <PerpSets,"PS">, <Elliptics,"El">, 
  <Hyperbolics,"Hy">];

for size := 2 to 4 do 
  SympSp := QuantumSymplecticSpace(size);

  for geometry in geometries do
    if not (geometry[1] eq Elliptics and size eq 2) then
      contextuality := [ IsContextual(instance[3]) : 
        instance in geometry[1](SympSp) ];
      
      notAllIdentical := contextuality[1];
      for value in contextuality do
        notAllIdentical xor:= value;
      end for;
      
      if notAllIdentical then
        Sprintf("%o(%o):INC(x%o)",geometry[2], size, #contextuality);
      else
        Sprintf("%o(%o):%o(x%o)",geometry[2], size, 
          contextuality[1] select "CON" else "NC ", #contextuality);
      end if;
    end if;
  end for;

end for;

exit;