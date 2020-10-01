> with(LinearAlgebra); with(SolveTools);

# DEFINES GENERAL FORM FOR 2x2x2x2x2 TENSORS
> form := sum(sum(sum(sum(sum(a[i1, i2, i3, i4, i5]*o[i1]*p[i2]*q[i3]*r[i4]*s[i5], i5 = 0 .. 1), i4 = 0 .. 1), i3 = 0 .. 1), i2 = 0 .. 1), i1 = 0 .. 1);

# CONVERT NUMBER IN BASIS 10 TO BINARY NUMBER
> dig2bin := proc (n, k) local N, L, i; N := n; L := [seq(0, i = 1 .. k)]; for i to k do L[i] := irem(N, 2); N := iquo(N, 2) end do; L end proc;

# FUNCTIONS USED FOR COMPUTING
> Formbis := proc (List) local S, i; S := {}; for i to 32 do S := {op(S), a[op(dig2bin(i-1, 5))] = List[i]} end do; S end proc;

> var := proc (aa, b, c, d, e) {o[`mod`(aa+1, 2)], p[`mod`(b+1, 2)], q[`mod`(c+1, 2)], r[`mod`(d+1, 2)], s[`mod`(e+1, 2)]} end proc;

> carte := proc (aa, b, c, d, e) local X; X := {o[aa] = 1, p[b] = 1, q[c] = 1, r[d] = 1, s[e] = 1} end proc;

> hessian_matrix := proc (F, aa, b, c, d, e) local A, tab, m, n; tab := [o[`mod`(aa+1, 2)], p[`mod`(b+1, 2)], q[`mod`(c+1, 2)], r[`mod`(d+1, 2)], s[`mod`(e+1, 2)]]; A := Matrix(5, 5); for m to 5 do for n to 5 do A(m, n) := diff(diff(F, tab[m]), tab[n]) end do end do; return A end proc;

> rank_hessian := proc (V, aa, b, c, d, e, val) local F; F := subs(carte(aa, b, c, d, e), subs(Formbis(V), form)); Rank(subs(val, hessian_matrix(F, aa, b, c, d, e))) end proc;


# FUNCTION THAT COMPUTES THE SINGULAR POINTS 
> Singbis := proc (V, par, cond) local singu, nsingu, eqns, aa, b, c, d, e, F, truc, ii; singu := {}; nsingu := {}; for aa from 0 to 1 do for b from 0 to 1 do for c from 0 to 1 do for d from 0 to 1 do for e from 0 to 1 do F := subs(carte(aa, b, c, d, e), subs(Formbis(V), form)); print(F); eqns := {op(cond), op({diff(F, o[0]) = 0, diff(F, o[1]) = 0, diff(F, p[0]) = 0, diff(F, p[1]) = 0, diff(F, q[0]) = 0, diff(F, q[1]) = 0, diff(F, r[0]) = 0, diff(F, r[1]) = 0, diff(F, s[0]) = 0, diff(F, s[1]) = 0, F = 0})}; truc := [solve(eqns, {op(var(aa, b, c, d, e)), op(par)})]; print(var(aa, b, c, d, e), [aa, b, c, d, e]); for ii to nops(truc) do select(has, truc[ii], {o[0], o[1], p[0], p[1], q[0], q[1], r[0], r[1], s[0], s[1]}); map(proc (x) options operator, arrow; op(2, x) end proc, %); if nops(select(has, %, {o[0], o[1], p[0], p[1], q[0], q[1], r[0], r[1], s[0], s[1]})) = 0 then singu := {op(singu), [eqns, truc[ii], rank_hessian(V, aa, b, c, d, e, truc[ii])]}; print(truc[ii], "Isolé", rank_hessian(V, a, b, c, d, e, truc[ii])) else nsingu := {op(nsingu), [eqns, truc[ii]]}; print(truc[ii], "Non isolé") end if end do; print(evalb(_SolutionsMayBeLost = true)) end do end do end do end do end do; save singu, `5_qubits_singu.txt`; save nsingu, `5_qubits_nsingu.txt` end proc;

# TESTING ON AN EXAMPLE
> Test := Normalize([1, 0, 0, sqrt(3), 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, sqrt(2), 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0], Euclidean);

> Singbis(Test, {}, {});


