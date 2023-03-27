/* S is PG(5,2)U{(0,...,0)} with the symplectic product named `InnerProduct` */
S := QSS(3);          
/* Sequence of tuples associated to hyperbolic quadrics in W(5,2) */
H := Hyperbolics(S);  
/* `points` is PG(5,2) */
points := { elt : elt in S | not IsZero(elt) }; 

hyperbolicHeptads := {};
for h in H do
  /* `hCompl` is the set of points that lie off the hyperbolic quadric `h` */
  hCompl := points diff h[2];  
  /* `elts` is `hCompl` ordered, in order to index each element */
  elts := {@ elt : elt in hCompl @}; 
  eltsIndex := {@ i : i in [1..#elts] @};
  edgesIndex := 
  { 
    { ij[1], ij[2] } 
      : ij in car<eltsIndex,eltsIndex> 
      | AntiCommute(elts[ij[1]],elts[ij[2]]) and ij[1] lt ij[2] 
  };
  graph := Graph<eltsIndex | edgesIndex>;
  heptads := 
  {
    { elts[Index(vertex)] : vertex in clc } 
      : clc in AllCliques(graph) 
      | #clc eq 7
  };  /* sets of seven observables that pairwise anticommute */
  hyperbolicHeptads join:= { <h,heptads> };
end for;

"Here are the heptads corresponding to a random hyperbolic of W(5,2), all the 
available heptads are stored in the variable `hyperbolicHeptads`:";
PrettyPrint(Random(hyperbolicHeptads)[2]);

heptads_ZZZ := {
  hHeptads[2]
    : hHeptads in hyperbolicHeptads 
    | hHeptads[1][1] eq S![1,0,1,0,1,0]
}; // the cardinal of this set is 1
computedQzzz := Random(heptads_ZZZ);
"Conwell heptads for Q+_(ZZZ)(5,2):";
PrettyPrint(computedQzzz);

computedQzzzString := {
 { ToPauliString(e) : e in h } : h in computedQzzz
};
QzzzMetod := {
  {"XII","ZXZ","ZIX","YYX","YZX","ZXY","YXI"},
  {"XII","ZZX","ZXI","YXY","YXZ","ZYX","YIX"},
  {"IXI","XZZ","IZX","YYX","ZYX","XZY","XYI"},
  {"IXI","ZZX","XZI","XYY","XYZ","YZX","IYX"},
  {"IIX","XZZ","IXZ","YXY","ZXY","XYZ","XIY"},
  {"IIX","ZXZ","XIZ","XYY","XZY","YXZ","IXY"},
  {"XXX","YXI","IYX","IZX","ZXI","XIZ","XIY"},
  {"XXX","YIX","IXY","IXZ","XZI","ZIX","YIX"}
};
"The computed Conwell heptads for Q+_(ZZZ)(5,2) as predicted by Metod:";
QzzzMetod;
"The computed Conwell heptads for Q+_(ZZZ)(5,2) are as predicted by Metod:";
QzzzMetod eq computedQzzzString;

QzzzMetodCorrected := {
  {"XII","ZXZ","ZIX","YYX","YZX","ZXY","YXI"},
  {"XII","ZZX","ZXI","YXY","YXZ","ZYX","YIX"},
  {"IXI","XZZ","IZX","YYX","ZYX","XZY","XYI"},
  {"IXI","ZZX","XZI","XYY","XYZ","YZX","IYX"},
  {"IIX","XZZ","IXZ","YXY","ZXY","XYZ","XIY"},
  {"IIX","ZXZ","XIZ","XYY","XZY","YXZ","IXY"},
  {"XXX","YXI","IYX","IZX","ZXI","XIZ","XIY"},
  {"XXX","XYI","IXY","IXZ","XZI","ZIX","YIX"}
};

"The computed Conwell heptads for Q+_(ZZZ)(5,2) are as predicted by Metod with
a correction (the first IYX is replaced with XYI):";
QzzzMetodCorrected eq computedQzzzString;

heptads_XXX := {
  hHeptads[2]
    : hHeptads in hyperbolicHeptads 
    | hHeptads[1][1] eq S![0,1,0,1,0,1]
}; // the cardinal of this set is 1
computedQxxx := Random(heptads_XXX);
"Conwell heptads for Q+_(XXX)(5,2):";
PrettyPrint(computedQxxx);

computedQxxxString := {
 { ToPauliString(e) : e in h } : h in computedQxxx
};
QxxxMetod := {
  {"ZII","XZX","XIZ","YYZ","YXZ","XZY","YZI"},
  {"ZII","XXZ","XZI","YZY","YZX","XYZ","YIZ"},
  {"IZI","ZXX","IXZ","YYZ","XYZ","ZXY","ZYI"},
  {"IZI","XXZ","ZXI","ZYY","ZYX","YXZ","IYZ"},
  {"IIZ","XZX","ZIX","ZYY","ZXY","YZX","IZY"},
  {"IIZ","ZXX","IZX","YZY","XZY","ZYX","ZIY"},
  {"ZZZ","YZI","IYZ","IXZ","XZI","ZIX","ZIY"},
  {"ZZZ","YIZ","IZY","IZX","ZXI","XIZ","YIZ"}
};
"The computed Conwell heptads for Q+_(XXX)(5,2) as predicted by Metod:";
QxxxMetod;
"The computed Conwell heptads for Q+_(XXX)(5,2) are as predicted by Metod:";
QxxxMetod eq computedQxxxString;

QxxxMetodCorrected := {
  {"ZII","XZX","XIZ","YYZ","YXZ","XZY","YZI"},
  {"ZII","XXZ","XZI","YZY","YZX","XYZ","YIZ"},
  {"IZI","ZXX","IXZ","YYZ","XYZ","ZXY","ZYI"},
  {"IZI","XXZ","ZXI","ZYY","ZYX","YXZ","IYZ"},
  {"IIZ","XZX","ZIX","ZYY","ZXY","YZX","IZY"},
  {"IIZ","ZXX","IZX","YZY","XZY","ZYX","ZIY"},
  {"ZZZ","YZI","IYZ","IXZ","XZI","ZIX","ZIY"},
  {"ZZZ","ZYI","IZY","IZX","ZXI","XIZ","YIZ"}
};

"The computed Conwell heptads for Q+_(XXX)(5,2) are as predicted by Metod with
a correction (the first YIZ is replaced with ZYI):";
QxxxMetodCorrected eq computedQxxxString;

// exit;
