restart;
with(LinearAlgebra);
with(VectorCalculus);
with(CodeGeneration);
with(FileTools);
TEST := false;
# State und Control Vektor
x := [r[1], r[2], r[3], q[1], q[2], q[3], q[4], v[1], v[2], v[3], omega[1], omega[2], omega[3], u[1], u[2], u[3], u[4]];

eq := Vector([q[1]^2+q[2]^2+q[3]^2+q[4]^2-1]);
eqMatrix := convert(eq, Matrix);

getNZeros := proc (M) local res, i, j, m, n; res := 0; i := 1; j := 1; m := 0; n := 0; m, n := Dimension(M); for i to m do for j to n do if M[i][j] <> 0 then res := res+1 end if end do end do; return res end proc;


# 

# Berechnung der nicht Nulleinträge der Hessematrizen
getNZerosHesses := proc (M, y) local m, n, CountHesse, k; CountHesse := 0; k := 1; m := RowDimension(M); for k to m do CountHesse := CountHesse+getNZeros(Hessian(M[k], y)) end do; return CountHesse end proc;


# Jacobimatrizen und Hessematrizen in Array umwandeln
getArray := proc (F, x) local m, n, JBtmp, CountJacobi, CountHesse, J, H, i, j, tmp, l, tmp2, h, k; JBtmp := Jacobian(F, x); CountJacobi := getNZeros(JBtmp); CountHesse := getNZerosHesses(F, x); J := Matrix(CountJacobi, 3); H := Matrix(CountHesse, 4); m, n := Dimension(JBtmp); k := 1; h := 1; for i to m do if F[i] <> 0 then for j to n do tmp := diff(F[i], x[j]); if tmp <> 0 then J[k, 1] := i; J[k, 2] := j; J[k, 3] := tmp; for l to n do tmp2 := diff(tmp, x[l]); if tmp2 <> 0 then H[h, 1] := i; H[h, 2] := j; H[h, 3] := l; H[h, 4] := tmp2; h := h+1 end if end do; k := k+1 end if end do end if end do; return J, H, CountJacobi, CountHesse end proc;

# Ausgabe Hesse und Jacobi in Array
J, H, CountJacobi, CountHesse := getArray(eq, x);

# Temporary Ordner finden
tmpDir := TemporaryDirectory();
# Temporary Ordner als aktuellen Ordner setzen.
currentdir(tmpDir);
currentdir();
# Optimieren in MATLAB und exportieren.
ConstraintsCount, m := Dimension(eqMatrix);

if TEST = false then Matlab(ConstraintsCount, resultname = "CountConstraints", defaulttype = integer, output = tmpRTConstraintsCount); Matlab(CountJacobi, resultname = "CountJacobi", defaulttype = integer, output = tmpRTConstraintsJacobiCount); Matlab(CountHesse, resultname = "CountHesse", defaulttype = integer, output = tmpRTConstraintsHesseCount); Matlab(eval(([codegen:-optimize])(eqMatrix, tryhard)), defaulttype = integer, output = tmpRTConstraintsFunction); Matlab(eval(([codegen:-optimize])(J, tryhard)), defaulttype = integer, output = tmpRTConstraintsJacobi)*Matlab(eval(([codegen:-optimize])(H, tryhard)), defaulttype = integer, output = tmpRTConstraintsHesse) else Matlab(eval(([codegen:-optimize])(eqMatrix, tryhard)), defaulttype = integer) end if;

# 
