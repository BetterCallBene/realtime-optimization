restart;
with(LinearAlgebra);
with(VectorCalculus);
with(CodeGeneration);
with(FileTools);
TEST := false;
DEB := false;

# State und Control Vektor
x := [r[1], r[2], r[3], q[1], q[2], q[3], q[4], v[1], v[2], v[3], omega[1], omega[2], omega[3], u[1], u[2], u[3], u[4]];

# Rotationsmatrix
R := proc (q) options operator, arrow; Matrix(3, 3, {(1, 1) = 1-2*q[3]^2-2*q[4]^2, (1, 2) = -2*q[1]*q[4]+2*q[2]*q[3], (1, 3) = 2*q[1]*q[3]+2*q[2]*q[4], (2, 1) = 2*q[1]*q[4]+2*q[2]*q[3], (2, 2) = 1-2*q[2]^2-2*q[4]^2, (2, 3) = -2*q[1]*q[2]+2*q[3]*q[4], (3, 1) = -2*q[1]*q[3]+2*q[2]*q[4], (3, 2) = 2*q[1]*q[2]+2*q[3]*q[4], (3, 3) = 1-2*q[2]^2-2*q[3]^2}) end proc;

# Inverse der Massenmatrix
Minv := Matrix(6, 6, {(1, 1) = 1/m, (1, 2) = 0, (1, 3) = 0, (1, 4) = 0, (1, 5) = 0, (1, 6) = 0, (2, 1) = 0, (2, 2) = 1/m, (2, 3) = 0, (2, 4) = 0, (2, 5) = 0, (2, 6) = 0, (3, 1) = 0, (3, 2) = 0, (3, 3) = 1/m, (3, 4) = 0, (3, 5) = 0, (3, 6) = 0, (4, 1) = 0, (4, 2) = 0, (4, 3) = 0, (4, 4) = 1/Iges[1], (4, 5) = 0, (4, 6) = 0, (5, 1) = 0, (5, 2) = 0, (5, 3) = 0, (5, 4) = 0, (5, 5) = 1/Iges[2], (5, 6) = 0, (6, 1) = 0, (6, 2) = 0, (6, 3) = 0, (6, 4) = 0, (6, 5) = 0, (6, 6) = 1/Iges[3]});

# Kreuzproduktmatrix [q_[1:3]x]:=Matrix(3, 3, {(1, 1) = 0, (1, 2) = -q[3], (1, 3) = q[2], (2, 1) = q[3], (2, 2) = 0, (2, 3) = -q[1], (3, 1) = -q[2], (3, 2) = q[1], (3, 3) = 0});
# 
getCrossProductMatrix := proc (q) local Omega; Omega := Matrix(3, 3, shape = antisymmetric); Omega[1, 2] := -q[3]; Omega[1, 3] := q[2]; Omega[2, 3] := -q[1]; return Omega end proc;

# q(x)p := [[[p[0],-p[1:3]^(T)],[p[1:3],p[0]I[3 x3]+[p[1:3]x]]]];
quatmultiply := proc (p, q) local pTmp, qTmp, M; pTmp := `<,>`(p[1], p[2], p[3], p[4]); qTmp := `<,>`(q[1], q[2], q[3], q[4]); M := Matrix(4); M[1] := pTmp[1]; M[1, 2 .. 4] := -Transpose(pTmp[2 .. 4]); M[2 .. 4, 1] := pTmp[2 .. 4]; M[2 .. 4, 2 .. 4] := pTmp[1]*IdentityMatrix(3)+getCrossProductMatrix(pTmp[2 .. 4]); return Multiply(M, qTmp) end proc;
# Theta(q, v, omega):=[[[m*omega x v + R(q)( )^(T) * [[[0,],[0,],[g,]]]],[omega x I[ges]omega]]];
Theta := proc (q, v, omega) local res; res := Vector(6); res[1 .. 3] := m*CrossProduct(omega, v)+m*Multiply(Transpose(R(q)), `<,>`(0, 0, g)); res[4 .. 6] := CrossProduct(omega, `<,>`(Iges[1]*omega[1], Iges[2]*omega[2], Iges[3]*omega[3])); return res end proc;
# T(omega, u):= [[[0,],[0,],[-k[T]*(&sum;)u[i]^(2),],[M*u^(2)+[[[0,],[0,],[1,]]],x omega *I[M]*(u[1]-u[2]+u[3]-u[4])]]];
T := proc (omega, u) local res, M; res := Vector(6); M := Matrix(3, 4, {(1, 1) = 0, (1, 2) = d*kT, (1, 3) = 0, (1, 4) = -d*kT, (2, 1) = -d*kT, (2, 2) = 0, (2, 3) = d*kT, (2, 4) = 0, (3, 1) = -kQ, (3, 2) = kQ, (3, 3) = -kQ, (3, 4) = kQ}); res[1 .. 3] := `<,>`(0, 0, kT*(u[1]^2+u[2]^2+u[3]^2+u[4]^2)); res[4 .. 6] := Multiply(M, `<,>`(u[1]^2, u[2]^2, u[3]^2, u[4]^2))+CrossProduct(`<,>`(0, 0, 1), omega)*IM*(u[1]-u[2]+u[3]-u[4]); return res end proc;

# Massenmatrix B(x)
B := proc (x) local vx, res; vx := Vector(x); res := Matrix(IdentityMatrix(13)); res[4 .. 7, 4 .. 7] := DiagonalMatrix(vx[4 .. 7]); return res end proc;
;


# dot(x):= [[[R(q)*v,],[1/(2)*q*[[[0,],[omega,]]],],[M^(-1)*(T(omega, u) - Theta(q, v, omega)),]]];
getDot := proc (x) local vx, res, tmpQ, tmpV, tmpOmega, tmpU; vx := Vector(x); res := Vector(13); tmpQ := vx[4 .. 7]; tmpV := vx[8 .. 10]; tmpOmega := vx[11 .. 13]; tmpU := vx[14 .. 17]; res[1 .. 3] := Multiply(R(tmpQ), tmpV); res[4 .. 7] := (1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega)); res[8 .. 13] := Multiply(Minv, T(tmpOmega, tmpU)-Theta(tmpQ, tmpV, tmpOmega)); return res end proc;
getDotTilde := proc (x) local vx, res, tmpQ, tmpV, tmpOmega, tmpU; vx := Vector(x); res := Vector(13); tmpQ := vx[4 .. 7]; tmpV := vx[8 .. 10]; tmpOmega := vx[11 .. 13]; tmpU := vx[14 .. 17]; res[1 .. 3] := Multiply(R(tmpQ), tmpV); res[4] := (-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[2]*q[2]-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[3]*q[3]-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[4]*q[4])[1]; res[5] := (-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[1]*q[1]-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[3]*q[3]-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[4]*q[4])[1]; res[6] := (-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[1]*q[1]-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[2]*q[2]-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[4]*q[4])[1]; res[7] := (-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[1]*q[1]-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[2]*q[2]-(1/2)*quatmultiply(tmpQ, `<,>`(0, tmpOmega))[3]*q[3])[1]; res[8 .. 13] := Multiply(Minv, T(tmpOmega, tmpU)-Theta(tmpQ, tmpV, tmpOmega)); return res end proc;

#  Ausführen von getDot(x)
dot := getDot(x);
dotTilde := getDotTilde(x);
dot := simplify(dot, 'symbolic');
dotTilde := simplify(dotTilde, 'symbolic');
# Für den Matlabexport in eine Matrix umwandeln

dotMatrix := convert(dot, Matrix);
dotMatrixTilde := convert(dotTilde, Matrix);
if DEB = true then dotMatrix[1 .. 3]; dotMatrix[4 .. 7]; dotMatrix[8 .. 13] end if;




if DEB = true then print(TestJ[11 .. 13, 14 .. 17]) end if;
# Berechnung der nicht Nulleinträge von Matrizen

getNZeros := proc (M) local res, i, j, m, n; res := 0; i := 1; j := 1; m := 0; n := 0; m, n := Dimension(M); for i to m do for j to n do if M[i][j] <> 0 then res := res+1 end if end do end do; return res end proc;


# 

# Berechnung der nicht Nulleinträge der Hessematrizen
getNZerosHesses := proc (M, y) local m, n, CountHesse, k; CountHesse := 0; k := 1; m := RowDimension(M); for k to m do CountHesse := CountHesse+getNZeros(Hessian(M[k], y)) end do; return CountHesse end proc;


# Jacobimatrizen und Hessematrizen in Array umwandeln
getArray := proc (F, x) local m, n, JBtmp, CountJacobi, CountHesse, J, H, i, j, tmp, l, tmp2, h, k; JBtmp := Jacobian(F, x); CountJacobi := getNZeros(JBtmp); CountHesse := getNZerosHesses(F, x); J := Matrix(CountJacobi, 3); H := Matrix(CountHesse, 4); m, n := Dimension(JBtmp); k := 1; h := 1; for i to m do if F[i] <> 0 then for j to n do tmp := diff(F[i], x[j]); if tmp <> 0 then J[k, 1] := i; J[k, 2] := j; J[k, 3] := tmp; for l to n do tmp2 := diff(tmp, x[l]); if tmp2 <> 0 then H[h, 1] := i; H[h, 2] := j; H[h, 3] := l; H[h, 4] := tmp2*ones; h := h+1 end if end do; k := k+1 end if end do end if end do; return J, H end proc;

# Ausgabe Hesse und Jacobi in Array
J, H := getArray(dot, x);


# 
JTilde, HTilde := getArray(dotTilde, x);
# Temporary Ordner finden
tmpDir := TemporaryDirectory();
# Temporary Ordner als aktuellen Ordner setzen.
currentdir(tmpDir);
currentdir();
# Optimieren in MATLAB und exportieren.

if TEST = false then Matlab(eval(([codegen:-optimize])(dotMatrix, tryhard)), defaulttype = integer, output = tmpRTOptFunction); VectorCalculus:-`*`(VectorCalculus:-`*`(VectorCalculus:-`*`(VectorCalculus:-`*`(Matlab(eval(([codegen:-optimize])(J, tryhard)), defaulttype = integer, output = tmpRTOptJacobi), Matlab(eval(([codegen:-optimize])(H, tryhard)), defaulttype = integer, output = tmpRTOptHesse)), Matlab(eval(([codegen:-optimize])(dotMatrixTilde, tryhard)), defaulttype = integer, output = tmpRTOptFTilde)), Matlab(eval(([codegen:-optimize])(JTilde, tryhard)), defaulttype = integer, output = tmpRTOptJTilde)), Matlab(eval(([codegen:-optimize])(HTilde, tryhard)), defaulttype = integer, output = tmpRTOptHTilde)) else if DEB = false then Matlab(dotMatrix, defaulttype = integer) else  end if end if;



