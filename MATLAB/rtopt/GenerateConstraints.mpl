restart;
with(LinearAlgebra);
with(VectorCalculus);
with(CodeGeneration);
with(FileTools);
TEST := false;
R := 3;
eps := 0.1e-1;
# State und Control Vektor
x := [r[1], r[2], r[3], q[1], q[2], q[3], q[4], v[1], v[2], v[3], omega[1], omega[2], omega[3], u[1], u[2], u[3], u[4]];

eqConstraints := Vector([VectorCalculus:-`+`(VectorCalculus:-`+`(VectorCalculus:-`+`(VectorCalculus:-`+`(q[1]^2, q[2]^2), q[3]^2), q[4]^2), -1)]);
inEqConstraints := Vector(`<,>`(VectorCalculus:-`+`(VectorCalculus:-`+`(R, VectorCalculus:-`-`(eps))^2, VectorCalculus:-`-`(VectorCalculus:-`+`(r[1]^2, r[2]^2))), VectorCalculus:-`+`(VectorCalculus:-`+`(r[1]^2, r[2]^2), VectorCalculus:-`-`(VectorCalculus:-`+`(R, eps)^2))));
eqMatrix := convert(eqConstraints, Matrix);

inEqMatrix := convert(inEqConstraints, Matrix);
getNZeros := proc (M) local res, i, j, m, n; res := 0; i := 1; j := 1; m := 0; n := 0; m, n := Dimension(M); for i to m do for j to n do if M[i][j] <> 0 then res := res+1 end if end do end do; return res end proc;


# 

# Berechnung der nicht Nulleinträge der Hessematrizen
getNZerosHesses := proc (M, y) local m, n, CountHesse, k; CountHesse := 0; k := 1; m := RowDimension(M); for k to m do CountHesse := CountHesse+getNZeros(Hessian(M[k], y)) end do; return CountHesse end proc;


# Jacobimatrizen und Hessematrizen in Array umwandeln
getArray := proc (F, x) local m, n, JBtmp, CountJacobi, CountHesse, J, H, i, j, tmp, l, tmp2, h, k; JBtmp := Jacobian(F, x); CountJacobi := getNZeros(JBtmp); CountHesse := getNZerosHesses(F, x); J := Matrix(CountJacobi, 3); H := Matrix(CountHesse, 4); m, n := Dimension(JBtmp); k := 1; h := 1; for i to m do if F[i] <> 0 then for j to n do tmp := diff(F[i], x[j]); if tmp <> 0 then J[k, 1] := i; J[k, 2] := j; J[k, 3] := tmp; for l to n do tmp2 := diff(tmp, x[l]); if tmp2 <> 0 then H[h, 1] := i; H[h, 2] := j; H[h, 3] := l; H[h, 4] := tmp2; h := h+1 end if end do; k := k+1 end if end do end if end do; return J, H, CountJacobi, CountHesse end proc;

# Ausgabe Hesse und Jacobi in Array
eqJ, eqH, eqCountJacobi, eqCountHesse := getArray(eqConstraints, x);
inEqJ, inEqH, inEqCountJacobi, inEqCountHesse := getArray(inEqConstraints, x);

# Temporary Ordner finden
tmpDir := TemporaryDirectory();
# Temporary Ordner als aktuellen Ordner setzen.
currentdir(tmpDir);
currentdir();
# Optimieren in MATLAB und exportieren.
eqCountConstraints, m := Dimension(eqMatrix);
inEqCountConstraints, m := Dimension(inEqMatrix);
if TEST = false then Matlab(eqCountConstraints, resultname = "eqCountConstraints", defaulttype = integer, output = tmpRTEqCountConstraints); Matlab(eqCountJacobi, resultname = "eqCountJacobi", defaulttype = integer, output = tmpRTEqCountConstraintsJacobi); Matlab(eqCountHesse, resultname = "eqCountHesse", defaulttype = integer, output = tmpRTEqCountConstraintsHesse); Matlab(eval(([codegen:-optimize])(eqMatrix, tryhard)), defaulttype = integer, output = tmpRTEqConstraintsFunction); VectorCalculus:-`*`(VectorCalculus:-`*`(Matlab(eval(([codegen:-optimize])(eqJ, tryhard)), defaulttype = integer, output = tmpRTEqConstraintsJacobi), Matlab(eval(([codegen:-optimize])(eqH, tryhard)), defaulttype = integer, output = tmpRTEqConstraintsHesse)), Matlab(inEqCountConstraints, resultname = "inEqCountConstraints", defaulttype = integer, output = tmpRTInEqCountConstraints)); Matlab(inEqCountJacobi, resultname = "inEqCountJacobi", defaulttype = integer, output = tmpRTInEqCountConstraintsJacobi); Matlab(inEqCountHesse, resultname = "inEqCountHesse", defaulttype = integer, output = tmpRTInEqCountConstraintsHesse); Matlab(eval(([codegen:-optimize])(inEqMatrix, tryhard)), defaulttype = integer, output = tmpRTInEqConstraintsFunction); VectorCalculus:-`*`(Matlab(eval(([codegen:-optimize])(inEqJ, tryhard)), defaulttype = integer, output = tmpRTInEqConstraintsJacobi), Matlab(eval(([codegen:-optimize])(inEqH, tryhard)), defaulttype = integer, output = tmpRTInEqConstraintsHesse)) else Matlab(eval(([codegen:-optimize])(eqMatrix, tryhard)), defaulttype = integer) end if;

# 
