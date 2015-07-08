function testCase()
A=randn(10,10); B=randn(10,10);
tic
C2 = tprod(A,[1 -2],B,[-2 2]);
toc
tic 
C  = etprod('ij',A,'ik',B,'kj'); mad(C2,C)
toc
C  = etprod({'i' 'j'},A,{'i' 'k'},B,{'k' 'j'});
A=randn(243,243, 13);B=randn(243,243,13);
tic
C3 = tprod(A,[-1 -2 1],B,[-1 -2 1])
toc
%C3 = tprod(B,[-1 -2 1],A,[-1 -2]);
%C3 = tprod(A,[-1 -2],B,[-1 -2 1]);
%C  = etprod('3',A,'12',B,'123');
%C  = etprod([3],A,[1 2],B,[1 2 3]);
%C  = etprod({'3'}, A,{'1' '2'},B,{'1' '2' '3'})
