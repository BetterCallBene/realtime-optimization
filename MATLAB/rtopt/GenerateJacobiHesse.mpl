restart;
with(LinearAlgebra);
with(VectorCalculus);
with(CodeGeneration);
with(FileTools);
x := [omega[1], omega[2], omega[3], u[1], u[2], u[3], u[4]];

T := proc (omega, u) options operator, arrow; (Vector(6, {(1) = 0, (2) = 0, (3) = kT*(sum(u[i]^2, i = 1 .. 4)), (4) = kT*d*(u[2]^2-u[4]^2), (5) = kT*d*(u[3]^2-u[1]^2), (6) = kQ*(-u[1]^2+u[2]^2-u[3]^2+u[4]^2)}))+(Vector(6, {(1) = 0, (2) = 0, (3) = 0, (4) = IM*omega[2]*(-u[1]+u[2]-u[3]+u[4]), (5) = IM*omega[1]*(u[1]-u[2]+u[3]-u[4]), (6) = 0})) end proc;
dynM := T(omega, u);


getNZeros := proc (M) local res, i, j, m, n; res := 0; i := 1; j := 1; m := 0; n := 0; m, n := Dimension(M); for i to m do for j to n do if M[i][j] <> 0 then res := res+1 end if end do end do; return res end proc;


getNZerosHesses := proc (M, y) local m, n, CountHesse, k; CountHesse := 0; k := 1; m := RowDimension(M); for k to m do CountHesse := CountHesse+getNZeros(Hessian(M[k], y)) end do; return CountHesse end proc;


getArray := proc (F, x) local m, n, JBtmp, CountJacobi, CountHesse, J, H, i, j, tmp, l, tmp2, h, k; JBtmp := Jacobian(F, x); CountJacobi := getNZeros(JBtmp); CountHesse := getNZerosHesses(F, x); J := Matrix(CountJacobi, 3); H := Matrix(CountHesse, 4); m, n := Dimension(JBtmp); k := 1; h := 1; for i to m do if F[i] <> 0 then for j to n do tmp := diff(F[i], x[j]); if tmp <> 0 then J[k, 1] := i; J[k, 2] := j; J[k, 3] := tmp; for l to n do tmp2 := diff(tmp, x[l]); if tmp2 <> 0 then H[h, 1] := i; H[h, 2] := j; H[h, 3] := l; H[h, 4] := tmp2; h := h+1 end if end do; k := k+1 end if end do end if end do; return J, H end proc;

J, H := getArray(dynM, x);
tmpDir := TemporaryDirectory();
currentdir(tmpDir);
Matlab(dynM, optimize, defaulttype = integer, output = tmpRTOptFunction);
Matlab(J, optimize, defaulttype = integer, output = tmpRTOptJacobi);
Matlab(H, optimize, defaulttype = integer, output = tmpRTOptHesse);

