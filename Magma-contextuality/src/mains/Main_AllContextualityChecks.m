geometries := [WLines, WBlocks, PerpSets, Elliptics, Hyperbolics];

for size := 2 to 4 do 
  SympSp := QuantumSymplecticSpace(size);

  for geometry in geometries do
    if not (geometry eq Elliptics and size eq 2) then
      geometryName := Split(Sprint(geometry,"Minimal"),"'")[2];
      contextuality := [ IsContextual(instance[3]) : instance in geometry(SympSp) ];
      notAllIdentical := contextuality[1];
      for value in contextuality do
        notAllIdentical xor:= value;
      end for;
      if notAllIdentical then
        Sprintf("%o on %o qubits are inconstant.",geometryName, size);
      else
        Sprintf("%o on %o qubits are all %o.",geometryName, size, 
          contextuality[1] select "contextual" else "non contextual");
      end if;
    end if;
  end for;

end for;

exit;