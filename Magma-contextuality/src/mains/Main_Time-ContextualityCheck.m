size := 3; 
SympSp := QuantumSymplecticSpace(size);

geometries := []
  cat WLines(SympSp) 
  cat WBlocks(SympSp) 
  cat PerpSets(SympSp) 
  cat ((size gt 2) select Elliptics(SympSp) else []) 
  cat Hyperbolics(SympSp)
  ;

// Time() seems to be broken, so we use Realtime()
t := Realtime();

contextuality := [ IsContextual(geometry[3]) : geometry in geometries ];

// for geometry in geometries do
//   c := IsContextual(geometry[3]);
// end for;

Realtime(t);
// all those contextuality checks take 4.4s for 3 qubits
// all those contextuality checks take 0.08s for 2 qubits

exit;