restart;
with(LinearAlgebra);
with(VectorCalculus);
with(CodeGeneration);
with(FileTools);
TEST := false;


x := [r[1], r[2], r[3], q[1], q[2], q[3], q[4], v[1], v[2], v[3], omega[1], omega[2], omega[3], u[1], u[2], u[3], u[4]];
xv := Vector(x);
R := proc (q) options operator, arrow; Matrix(3, 3, {(1, 1) = 1-2*q[3]^2-2*q[4]^2, (1, 2) = -2*q[1]*q[4]+2*q[2]*q[3], (1, 3) = 2*q[1]*q[3]+2*q[2]*q[4], (2, 1) = 2*q[1]*q[4]+2*q[2]*q[3], (2, 2) = 1-2*q[2]^2-2*q[4]^2, (2, 3) = -2*q[1]*q[2]+2*q[3]*q[4], (3, 1) = -2*q[1]*q[3]+2*q[2]*q[4], (3, 2) = 2*q[1]*q[2]+2*q[3]*q[4], (3, 3) = 1-2*q[2]^2-2*q[3]^2}) end proc;
R3 := proc (q) options operator, arrow; Vector(3, {(1) = -2*q[1]*q[3]+2*q[2]*q[4], (2) = 2*q[1]*q[2]+2*q[3]*q[4], (3) = 1-2*q[2]^2-2*q[3]^2}) end proc;
Minv := Matrix(6, 6, {(1, 1) = 1/m, (1, 2) = 0, (1, 3) = 0, (1, 4) = 0, (1, 5) = 0, (1, 6) = 0, (2, 1) = 0, (2, 2) = 1/m, (2, 3) = 0, (2, 4) = 0, (2, 5) = 0, (2, 6) = 0, (3, 1) = 0, (3, 2) = 0, (3, 3) = 1/m, (3, 4) = 0, (3, 5) = 0, (3, 6) = 0, (4, 1) = 0, (4, 2) = 0, (4, 3) = 0, (4, 4) = 1/Iges[1], (4, 5) = 0, (4, 6) = 0, (5, 1) = 0, (5, 2) = 0, (5, 3) = 0, (5, 4) = 0, (5, 5) = 1/Iges[2], (5, 6) = 0, (6, 1) = 0, (6, 2) = 0, (6, 3) = 0, (6, 4) = 0, (6, 5) = 0, (6, 6) = 1/Iges[3]});
Q := Typesetting[delayDotProduct](1/2, Vector(4, {(1) = -q[2]*omega[1]-q[3]*omega[2]-q[4]*omega[3], (2) = q[1]*omega[1]+q[3]*omega[3]-q[4]*omega[2], (3) = q[1]*omega[2]-q[2]*omega[3]+q[4]*omega[1], (4) = q[1]*omega[3]+q[2]*omega[2]-q[3]*omega[1]}), true);
Theta := proc (q, v, omega) options operator, arrow; Vector(6, {(1) = m*(omega[2]*v[3]-omega[3]*v[2]+R3(q)[1]*g), (2) = m*(omega[3]*v[1]-omega[1]*v[3]+R3(q)[2]*g), (3) = m*(omega[1]*v[2]-omega[2]*v[1]+R3(q)[3]*g), (4) = Iges[3]*omega[2]*omega[3]-Iges[2]*omega[3]*omega[2], (5) = Iges[1]*omega[3]*omega[1]-Iges[3]*omega[1]*omega[3], (6) = Iges[2]*omega[1]*omega[2]-Iges[1]*omega[2]*omega[1]}) end proc;

T := proc (omega, u) options operator, arrow; (Vector(6, {(1) = 0, (2) = 0, (3) = kT*(sum(u[i]^2, i = 1 .. 4)), (4) = kT*d*(u[2]^2-u[4]^2), (5) = kT*d*(u[3]^2-u[1]^2), (6) = kQ*(-u[1]^2+u[2]^2-u[3]^2+u[4]^2)}))+(Vector(6, {(1) = 0, (2) = 0, (3) = 0, (4) = IM*omega[2]*(-u[1]+u[2]-u[3]+u[4]), (5) = IM*omega[1]*(u[1]-u[2]+u[3]-u[4]), (6) = 0})) end proc;
dot := Vector([Multiply(R(q), xv[8 .. 10]), Q, Multiply(Minv, VectorCalculus:-`+`(T(omega, u), VectorCalculus:-`-`(Theta(q, v, omega))))]);


dynM := dot;


getNZeros := proc (M) local res, i, j, m, n; res := 0; i := 1; j := 1; m := 0; n := 0; m, n := Dimension(M); for i to m do for j to n do if M[i][j] <> 0 then res := res+1 end if end do end do; return res end proc;


getNZerosHesses := proc (M, y) local m, n, CountHesse, k; CountHesse := 0; k := 1; m := RowDimension(M); for k to m do CountHesse := CountHesse+getNZeros(Hessian(M[k], y)) end do; return CountHesse end proc;


getArray := proc (F, x) local m, n, JBtmp, CountJacobi, CountHesse, J, H, i, j, tmp, l, tmp2, h, k; JBtmp := Jacobian(F, x); CountJacobi := getNZeros(JBtmp); CountHesse := getNZerosHesses(F, x); J := Matrix(CountJacobi, 3); H := Matrix(CountHesse, 4); m, n := Dimension(JBtmp); k := 1; h := 1; for i to m do if F[i] <> 0 then for j to n do tmp := diff(F[i], x[j]); if tmp <> 0 then J[k, 1] := i; J[k, 2] := j; J[k, 3] := tmp; for l to n do tmp2 := diff(tmp, x[l]); if tmp2 <> 0 then H[h, 1] := i; H[h, 2] := j; H[h, 3] := l; H[h, 4] := tmp2*ones; h := h+1 end if end do; k := k+1 end if end do end if end do; return J, H end proc;

J, H := getArray(dynM, x);

tmpDir := TemporaryDirectory();
currentdir(tmpDir);
currentdir();
if TEST = false then Matlab(dynM, optimize, defaulttype = integer, output = tmpRTOptFunction); VectorCalculus:-`*`(Matlab(eval(([codegen:-optimize])(J, tryhard)), defaulttype = integer, output = tmpRTOptJacobi), Matlab(eval(([codegen:-optimize])(H, tryhard)), defaulttype = integer, output = tmpRTOptHesse)) else Matlab(dynM, optimize, defaulttype = integer); Matlab(eval(([codegen:-optimize])(J, tryhard)), defaulttype = integer); Matlab(eval(([codegen:-optimize])(H, tryhard)), defaulttype = integer) end if;


