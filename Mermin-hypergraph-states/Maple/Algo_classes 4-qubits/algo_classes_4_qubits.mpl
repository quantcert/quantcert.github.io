> with(LinearAlgebra);



# DEFINITION OF THE INVARIANTS
> Operateur_A := proc (Y) local A; A := add(add(add(add(Y(i+2*j+4*k+8*l+1)*x[i+1]*y[j+1]*z[k+1]*t[l+1], i = 0 .. 1), j = 0 .. 1), k = 0 .. 1), l = 0 .. 1); return A end proc;

> Operateur_H := proc (Y) local n; n := Y(1)*Y(16)-Y(15)*Y(2)-Y(3)*Y(14)+Y(13)*Y(4)-Y(5)*Y(12)+Y(11)*Y(6)+Y(7)*Y(10)-Y(9)*Y(8); return n end proc;

> Operateur_L := proc (Y) local n, L; L := Matrix(4, 4, [[Y(1), Y(5), Y(9), Y(13)], [Y(2), Y(6), Y(10), Y(14)], [Y(3), Y(7), Y(11), Y(15)], [Y(4), Y(8), Y(12), Y(16)]]); n := Determinant(L); return n end proc;

> Operateur_M := proc (Y) local n, M; M := Matrix(4, 4, [[Y(1), Y(9), Y(3), Y(11)], [Y(2), Y(10), Y(4), Y(12)], [Y(5), Y(13), Y(7), Y(15)], [Y(6), Y(14), Y(8), Y(16)]]); n := Determinant(M); return n end proc;

> Operateur_N := proc (Y) local n; n := -Operateur_L(Y)-Operateur_M(Y); return n end proc;

> Operateur_Dxy := proc (Y) local A, b_xy, Dxy, B_xy; A := add(add(add(add(Y(i+2*j+4*k+8*l+1)*x[i+1]*y[j+1]*z[k+1]*t[l+1], i = 0 .. 1), j = 0 .. 1), k = 0 .. 1), l = 0 .. 1); b_xy := Determinant(Matrix(2, 2, [seq(seq(diff(diff(A, z[i]), t[j]), i = 1 .. 2), j = 1 .. 2)])); B_xy := Matrix(3, 3, [seq(seq(seq(seq(mvcoeff(b_xy, x[i]*x[j]*y[k]*y[l]), k = l .. 2), l = 1 .. 2), i = j .. 2), j = 1 .. 2)]); Dxy := -Determinant(B_xy); return Dxy end proc;

> Operateur_S := proc (Y) local n; n := (1/12)*Operateur_H(Y)^4-(2/3)*Operateur_H(Y)^2*Operateur_L(Y)+(2/3)*Operateur_H(Y)^2*Operateur_M(Y)-2*(Operateur_H(Y))(Operateur_Dxy(Y))+(4/3)*Operateur_L(Y)^2+(4/3)*Operateur_L(Y)*Operateur_M(Y)+(4/3)*Operateur_M(Y)^2; return n end proc;

> Operateur_T := proc (Y) local n; n := (1/216)*Operateur_H(Y)^6-(1/18)*Operateur_H(Y)^4*(Operateur_L(Y)-Operateur_M(Y))-(1/6)*Operateur_H(Y)^3*Operateur_Dxy(Y)+(1/9)*Operateur_H(Y)^2*(2*Operateur_L(Y)^2-Operateur_L(Y)*Operateur_M(Y)+2*Operateur_M(Y)^2)+(2/3)*Operateur_H(Y)*(Operateur_L(Y)-Operateur_M(Y))*Operateur_Dxy(Y)-(8/27)*Operateur_L(Y)^3+(8/27)*Operateur_M(Y)^3-(4/9)*Operateur_L(Y)*Operateur_M(Y)*(Operateur_L(Y)-Operateur_M(Y))+Operateur_Dxy(Y)^2; return n end proc;

> mvcoeff := proc (expr, term) local a, t; a := expand(expr); if type(term, `^`) or type(term, name) then return coeff(a, term) end if; for t in [op(term)] do a := coeff(a, t) end do; return a end proc;

> dig2bin := proc (n, k) local N, L, i; N := n; L := [seq(0, i = 1 .. k)]; for i to k do L[i] := irem(N, 2); N := iquo(N, 2) end do; L end proc;

> nonnul := proc (P) if expand(P) <> 0 then 1 else 0 end if end proc;

> VecEvalCov2 := proc (S, V) local V2; V2 := V; V22 := subs(S, V2); V3 := map(expand, V22); V4 := map(simplify, V3); V5 := map(nonnul, V4); return V5 end proc;

> Substitue := proc (Y) local i, S; S := {}; for i to 16 do S := {op(S), a[op(dig2bin(i-1, 4))] = Y(i)} end do end proc;

> IsNilpotent := proc (Y) local A; A := [Operateur_H(Y), Operateur_L(Y), Operateur_M(Y), Operateur_Dxy(Y)]; A := map(expand, A); A := map(evalf, A); evalb(A = [0., 0., 0., 0.]) end proc;

> Quartic_1 := proc (Y) local n; n := x^4-2*Operateur_H(Y)*x^3*y+(Operateur_H(Y)^2+2*Operateur_L(Y)+4*Operateur_M(Y))*x^2*y^2+(4*Operateur_Dxy(Y)-4*Operateur_H(Y)*(Operateur_M(Y)+(1/2)*Operateur_L(Y)))*x*y^3+Operateur_L(Y)^2*y^4; return n end proc;

> Quartic_2 := proc (Y) local n; n := x^4-2*Operateur_H(Y)*x^3*y+(Operateur_H(Y)^2-4*Operateur_L(Y)-2*Operateur_M(Y))*x^2*y^2+(-2*Operateur_M(Y)*Operateur_H(Y)+4*Operateur_Dxy(Y))*x*y^3+Operateur_M(Y)^2*y^4; return n end proc;

> Quartic_3 := proc (Y) local n; n := x^4-2*Operateur_H(Y)*x^3*y+(Operateur_H(Y)^2+2*Operateur_L(Y)-2*Operateur_M(Y))*x^2*y^2-(2*Operateur_L(Y)*Operateur_H(Y)+2*Operateur_M(Y)*Operateur_H(Y)-4*Operateur_Dxy(Y))*x*y^3+Operateur_N(Y)^2*y^4; return n end proc;

> Operateur_I2 := proc (Y) local n, L, B, M, Dxy; L := Operateur_L(Y); B := Operateur_H(Y); M := Operateur_M(Y); Dxy := Operateur_Dxy(Y); n := (4/3)*L^2-(4/3)*B^2*M+(4/3)*L*M+(4/3)*M^2+2*B*Dxy+(1/12)*B^4-(2/3)*B^2*L; return n end proc;

> Operateur_I3 := proc (Y) local n, L, B, M, Dxy; L := Operateur_L(Y); B := Operateur_H(Y); M := Operateur_M(Y); Dxy := Operateur_Dxy(Y); n := (8/27)*L^3-(1/216)*B^6-(8/27)*M^3-(1/6)*Dxy*B^3+(4/3)*B*M*Dxy-(5/9)*B^2*M*L+(2/3)*B*L*Dxy-(2/9)*B^2*L^2-(5/9)*B^2*M^2-Dxy^2+(4/9)*L^2*M+(1/18)*B^4*L+(1/9)*B^4*M-(4/9)*L*M^2; return n end proc;

> Hyperdeterminant := proc (Y) local n, H, L, M, Dxy, Q1; H := Operateur_H(Y); L := Operateur_L(Y); Dxy := Operateur_Dxy(Y); M := Operateur_M(Y); Q1 := x^4-2*H*x^3+(H^2+2*L+4*M)*x^2-(-4*Dxy+4*H*(M+(1/2)*L))*x+L^2; n := discrim(Q1, x); return n end proc;

> Hessian := proc (f) local hess; hess := Determinant(Matrix(2, 2, [diff(diff(f, x), x), diff(diff(f, x), y), diff(diff(f, y), x), diff(diff(f, y), y)])); return hess end proc;

> Jacobian_of_Hessian := proc (f) local hess, jacob; hess := Hessian(f); jacob := Determinant(Matrix(2, 2, [diff(f, x), diff(f, y), diff(hess, x), diff(hess, y)])); return jacob end proc;



##COVARIANTS DEFINITION

# COVARIANTS USED FOR NON-NILPOTENT FORMS
> read "vect_cov1.txt";

# COVARIANTS USED FOR NILPOTENT FORMS
> read "vect_cov2.txt";


# DETERMINE THE NILPOTENT STRATA
> Nilpotent_type := proc (Y) local A, Subs, Eval; A := Operateur_A(Y); if IsNilpotent(Y) then if A = 0 then return [Gr[0]] else Subs := Substitue(Y); Eval := VecEvalCov2(Subs, Vect_cov_nilpotent_cone); if Eval = [0, 0, 0, 0, 0, 0, 0] then return [Gr[1]] elif Eval = [1, 0, 0, 0, 0, 0, 0] then return [Gr[2]] elif Eval = [1, 1, 0, 0, 0, 0, 0] then return [Gr[3]] elif Eval = [1, 1, 0, 1, 0, 0, 0] then return [Gr[4]] elif Eval = [1, 1, 1, 0, 0, 0, 0] then return [Gr[5]] elif Eval = [1, 1, 1, 1, 1, 0, 0] then return [Gr[6]] elif Eval = [1, 1, 1, 1, 1, 1, 0] then return [Gr[7]] elif Eval = [1, 1, 1, 1, 1, 1, 1] then return [Gr[8]] end if end if else return "Not a Nilpotent form" end if end proc;

# DETERMINE THE VERSTRAETE FAMILY
> Verstraete_type := proc (Y) local hyper, Q1, Q2, Q3, T1, T2, T3, I2, I3, Hess1, Hess2, Hess3, LL, B, M, N, Dxy, S, V, Eval; L := 'L'; hyper := expand(Hyperdeterminant(Y)); Q1 := expand(Quartic_1(Y)); Q2 := expand(Quartic_2(Y)); Q3 := expand(Quartic_3(Y)); Hess1 := expand(Hessian(Q1)); Hess2 := expand(Hessian(Q2)); Hess3 := expand(Hessian(Q3)); T1 := expand(Jacobian_of_Hessian(Q1)); T2 := expand(Jacobian_of_Hessian(Q2)); T3 := expand(Jacobian_of_Hessian(Q3)); I2 := expand(Operateur_I2(Y)); I3 := expand(Operateur_I3(Y)); LL := expand(Operateur_L(Y)); B := expand(Operateur_H(Y)); M := expand(Operateur_M(Y)); N := expand(Operateur_N(Y)); Dxy := expand(Operateur_Dxy(Y)); if IsNilpotent(Y) then return Nilpotent_type(Y) end if; if LL = 0 and M = 0 then if Dxy = 0 and B = 0 then return "Nilpotent form" elif Dxy <> 0 and hyper <> 0 then return [G[abcd], d = 0] elif Dxy <> 0 and hyper = 0 then S := Substitue(Y); V := [Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0] then return [G[abcd], b = a, c = -2*a, d = 0] else return [L[abc[2]], a = 0, c = (1/2)*b] end if elif Dxy = 0 and B <> 0 then S := Substitue(Y); V := [Cov_Gbar, Cov_G, Cov_H, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0, 0, 0] then return [G[abcd], a = 0, b = 0, c = d] elif Eval = [0, 1, 1, 0] then return [L[abc[2]], a = b, c = 0] elif Eval = [0, 0, 1, 0] then return [L[abc[2]], a = 0, b = 0] elif Eval = [1, 1, 1, 0] then return [L[a[2]*b[2]], a = 0] elif Eval = [1, 1, 1, 1] then return [L[a[2]*O[`&oplus;`(3, conjugate(1))]]] else return "Error, unexpected case 1" end if end if elif LL = 0 and M <> 0 then if hyper <> 0 then return [G[abcd], d = 0] elif Dxy = B*M and B*B+4*M <> 0 then S := Substitue(Y); V := [Cov_K3, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0] then return [G[abcd], c = d, d = 0] elif Eval = [1, 0] then return [L[abc[2]], c = 0] elif Eval = [1, 1] then return [L[a[2]*b[2]]] else return "Error, unexpected case 2" end if elif Dxy <> B*M and T1 <> 0 and T2 <> 0 then S := Substitue(Y); V := [Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0] then return [G[abcd], b = c, d = 0] else return [L[abc[2]], b = 0] end if elif Dxy = B*M and B*B+4*M = 0 and Hess2 = 0 then S := Substitue(Y); V := [Cov_C, Cov_D, Cov_K5, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0, 0, 0] then return [G[abcd], b = 0, c = 0, d = 0] elif Eval = [1, 0, 0, 0] then return [L[abc[2]], b = 0, c = 0] elif Eval = [1, 1, 1, 0] then return [L[ab[3]], a = 0] elif Eval = [1, 1, 0, 0] then return [L[a[2]*b[2]], a = b] elif Eval = [1, 1, 1, 1] then return [L[a[4]]] else return "Error, unexpected case 3" end if elif Dxy <> B*M and I2 = 0 and hyper = 0 then if B = 0 then return "Nilpotent form" else S := Substitue(Y); V := [Cov_D, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0] then return [G[abcd], a = b, b = c, d = 0] elif Eval = [1, 0] then return [L[abcc], b = c, a = 0] elif Eval = [1, 1] then return [L[abbb], b = 0] else return "Error, unexpected case 4" end if end if else return "Error, unexpected case 5" end if elif M = 0 and LL <> 0 then if hyper <> 0 then return [G[abcd], d = 0] elif Dxy = 0 and B*B <> 4*LL then S := Substitue(Y); V := [Cov_K3, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0] then return [G[abcd], c = d, d = 0] elif Eval = [1, 0] then return [L[abc[2]], c = 0] elif Eval = [1, 1] then return [L[a[2]*b[2]]] else return "Error, unexpected case 6" end if elif Dxy <> 0 and T1 <> 0 and T2 <> 0 then S := Substitue(Y); V := [Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0] then return [G[abcd], b = c, d = 0] else return [L[abc[2]], b = 0] end if elif Dxy = 0 and B*B = 4*LL then S := Substitue(Y); V := [Cov_C, Cov_D, Cov_K5, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0, 0, 0] then return [G[abcd], b = 0, c = 0, d = 0] elif Eval = [1, 0, 0, 0] then return [L[abc[2]], b = 0, c = 0] elif Eval = [1, 1, 1, 0] then return [L[ab[3]], a = 0] elif Eval = [1, 1, 0, 0] then return [L[a[2]*b[2]], a = b] elif Eval = [1, 1, 1, 1] then return [L[a[4]]] else return "Error, unexpected case 7" end if elif Dxy <> 0 and I2 = 0 and hyper = 0 then if B = 0 then return "Nilpotent form" else S := Substitue(Y); V := [Cov_D, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0] then return [G[abcd], a = b, b = c, d = 0] elif Eval = [1, 0] then return [L[abcc], b = c, a = 0] elif Eval = [1, 1] then return [L[abbb], b = 0] else return "Error, unexpected case 8" end if end if else return "Error, unexpected case 9" end if elif N = 0 and LL <> 0 and M <> 0 then if hyper <> 0 then return [G[abcd], d = 0] elif Dxy = 0 and B*B <> 4*M then S := Substitue(Y); V := [Cov_K3, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0] then return [G[abcd], c = d, d = 0] elif Eval = [1, 0] then return [L[abc[2]], c = 0] elif Eval = [1, 1] then return [L[a[2]*b[2]]] else return "Error, unexpected case 10" end if elif Dxy <> 0 and T1 <> 0 and T2 <> 0 then S := Substitue(Y); V := [Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0] then return [G[abcd], b = c, d = 0] else return [L[abc[2]], b = 0] end if elif Dxy = 0 and B*B = 4*M then S := Substitue(Y); V := [Cov_C, Cov_D, Cov_K5, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0, 0, 0] then return [G[abcd], b = 0, c = 0, d = 0] elif Eval = [1, 0, 0, 0] then return [L[abc[2]], b = 0, c = 0] elif Eval = [1, 1, 1, 0] then return [L[ab[3]], a = 0] elif Eval = [1, 1, 0, 0] then return [L[a[2]*b[2]], a = b] elif Eval = [1, 1, 1, 1] then return [L[a[4]]] else return "Error, unexpected case 11" end if elif Dxy <> 0 and I2 = 0 and hyper = 0 then if B = 0 then return "Nilpotent form" else S := Substitue(Y); V := [Cov_D, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0] then return [G[abcd], a = b, b = c, d = 0] elif Eval = [1, 0] then return [L[abcc], b = c, a = 0] elif Eval = [1, 1] then return [L[abbb], b = 0] else return "Error, unexpected case 12" end if end if else return "Error, unexpected case 13" end if else if hyper <> 0 then return [G[abcd]] elif T1 <> 0 and T2 <> 0 and T3 <> 0 and I2 <> 0 and I3 <> 0 then S := Substitue(Y); V := [Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0] then return [G[abcd], c = d] else return [L[abcc]] end if elif I2 = 0 and I3 = 0 and Hess1 <> 0 then S := Substitue(Y); V := [Cov_K5, Cov_L]; Eval := VecEvalCov2(S, V); if Eval = [0, 0] then return [G[abcd], b = c, c = d] elif Eval = [1, 0] then return [L[abc[2]], b = c] elif Eval = [1, 1] then return [L[ab[3]]] else return "Error, unexpected case 14" end if else return "Error, unexpected case 15" end if end if end proc;



# EVALUATE THE CLASSIFICATION FUNCTIONS ON ALL NORMAL FORMS

> aaa := 'aaa'; bbb := 'bbb'; ccc := 'ccc'; ddd := 'ddd';

> Y1 := `<,>`((aaa+ddd)*(1/2), 0, 0, (aaa-ddd)*(1/2), 0, (bbb+ccc)*(1/2), (bbb-ccc)*(1/2), 0, 0, (bbb-ccc)*(1/2), (bbb+ccc)*(1/2), 0, (aaa-ddd)*(1/2), 0, 0, (aaa+ddd)*(1/2)); Verstraete_type(Y1); IsNilpotent(Y1); factor(Hyperdeterminant(Y1));

> Y2 := `<,>`((aaa+bbb)*(1/2), 0, 0, (aaa-bbb)*(1/2), 0, ccc, 1, 0, 0, 0, ccc, 0, (aaa-bbb)*(1/2), 0, 0, (aaa+bbb)*(1/2)); Verstraete_type(Y2); IsNilpotent(Y2); Hyperdeterminant(Y2);

> Y3 := `<,>`(aaa, 0, 0, 1, 0, bbb, 1, 0, 0, 0, bbb, 0, 0, 0, 0, aaa); Verstraete_type(Y3); IsNilpotent(Y3); Hyperdeterminant(Y3);

> Y4 := `<,>`(aaa, I/sqrt(2), I/sqrt(2), 0, 0, (aaa+bbb)*(1/2), (aaa-bbb)*(1/2), (-I)*(1/sqrt(2)), 0, (aaa-bbb)*(1/2), (aaa+bbb)*(1/2), (-I)*(1/sqrt(2)), 0, 0, 0, aaa); Verstraete_type(Y4); IsNilpotent(Y4); Hyperdeterminant(Y4);

> Y5 := `<,>`(aaa, I, 0, 0, 0, aaa, 1, 0, 0, 0, aaa, -I, 0, 0, 0, aaa); Verstraete_type(Y5); IsNilpotent(Y5); Hyperdeterminant(Y5);

> Y6 := `<,>`(aaa, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, aaa); Verstraete_type(Y6); IsNilpotent(Y6); Hyperdeterminant(Y6);

> Y7 := `<,>`(1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0); Verstraete_type(Y7); IsNilpotent(Y7); Hyperdeterminant(Y7);

> Y8 := `<,>`(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0); Verstraete_type(Y8); IsNilpotent(Y8); Hyperdeterminant(Y8);

> Y9 := `<,>`(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0); Verstraete_type(Y9); IsNilpotent(Y9); Hyperdeterminant(Y9);


> N0 := `<,>`(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

> N1 := `<,>`(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

> N2 := `<,>`(1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

> N3 := `<,>`(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0);

> N4 := `<,>`(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0);

> N5 := `<,>`(0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0);

> N6 := `<,>`(1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0);

> N7 := `<,>`(0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0);

> N8 := `<,>`(0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0);

> Nilpotent_type(N0); Nilpotent_type(N1); Nilpotent_type(N2); Nilpotent_type(N3); Nilpotent_type(N4); Nilpotent_type(N5); Nilpotent_type(N6); Nilpotent_type(N7); Nilpotent_type(N8);



# IMPLEMENTATION OF GROVER'S ALGORITHM

> Grover_algorithm_qubits := proc (Liste, nb_qubits, nb_iterations) local Y, i, j, NB, m, n, moyenne; m := nb_iterations; n := nb_qubits; Y := Vector(2^n, 1/sqrt(2^n)); NB := nops(Liste); for i to m do for j to NB do Y(Liste[j]) := -Y(Liste[j]) end do; moyenne := 0; for j to 2^n do moyenne := moyenne+Y(j) end do; moyenne := moyenne/2^n; for j to 2^n do Y(j) := 2*moyenne-Y(j) end do end do; return Y end proc;


> Grover_algorithm_qubits([2], 3, 1);

# PLOTTING THE HYPERDETERMINANT
> SeqX := []; SeqY := []; for GGG from 0 to 50 do Y := Grover_algorithm_qubits([1, 2, 3, 6, 11, 16], 4, GGG); SeqX := [op(SeqX), GGG]; SeqY := [op(SeqY), abs(Hyperdeterminant(Y))] end do; plot(SeqX, SeqY);


# DEFINITION OF GRAPH AND HYPERGRAPH STATES FOR 4-QUBITS
> G1 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1], Euclidean); G2 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, -1], Euclidean); G3 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, 1], Euclidean); G4 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1], Euclidean); G5 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1], Euclidean); G6 := Normalize([1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, -1], Euclidean); G7 := Normalize([1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, 1, -1, 1], Euclidean); G8 := Normalize([1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1], Euclidean); G9 := Normalize([1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1], Euclidean); G10 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1], Euclidean); G11 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1], Euclidean); G12 := Normalize([1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1], Euclidean); G13 := Normalize([1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, -1, -1], Euclidean); G14 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1], Euclidean); G15 := Normalize([1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1], Euclidean); G16 := Normalize([1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1], Euclidean); G17 := Normalize([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1], Euclidean); G18 := Normalize([1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1], Euclidean); G19 := Normalize([1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, -1], Euclidean); G20 := Normalize([1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, 1, 1, -1, 1, -1], Euclidean); G21 := Normalize([1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, 1], Euclidean); G22 := Normalize([1, 1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1], Euclidean); G23 := Normalize([1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1], Euclidean); G24 := Normalize([1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1], Euclidean); G25 := Normalize([1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, -1, 1, -1, 1, -1], Euclidean); G26 := Normalize([1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1], Euclidean); G27 := Normalize([1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1], Euclidean); G4_5 := Normalize([1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1], Euclidean); G4_4 := Normalize([1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1], Euclidean);


> Verstraete_type(G1); Hyperdeterminant(G1);


