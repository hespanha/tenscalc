% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

clear all

profile clear;
profile on;

code=csparse();

Tvariable a0 []
Tvariable b0 []
Tvariable a1 [3]
Tvariable b1 [3]
Tvariable c1 [7]
Tvariable a2 [7,3]
Tvariable b2 [7,7]
Tvariable a3 [3,7,11]
Tvariable b3 [3,7,11]

declareSet(code,a0,'a0');
declareSet(code,b0,'b0');
declareSet(code,a1,'a1');
declareSet(code,b1,'b1');
declareSet(code,a2,'a2');
declareSet(code,b2,'b2');
declareSet(code,a3,'a3');
declareSet(code,b3,'b3');

disp('plus, minus, reshape, subsref, ones, repmat, clp')
y0=a0+b0
yy0=-y0+1
y1=a1-b1
yy1=a1-5
y2=-a3+b3
yy2=4-b3
yy3=0-b3
yy4=repmat(a2,[2,3])
y5=clp(a0,b0)
y6=clp(a1,b1)
y7=clp(a3,b3)

declareGet(code,y0,'y0')
declareGet(code,yy0,'yy0')
declareGet(code,y1,'y1')
declareGet(code,yy1,'yy1')
declareGet(code,y2,'y2')
declareGet(code,yy2,'yy2')
declareGet(code,yy3,'yy3')
declareGet(code,yy4,'yy4')
declareGet(code,y5,'y5')
declareGet(code,y6,'y6')
declareGet(code,y7,'y7')

disp('==, <=, >= , <, >=')
q01=a0==b0
q02=a0>=b0
q03=a0<=b0
q04=a0>b0
q05=a0<b0
q11=a1==b1
q12=a1>=b1
q13=a1<=b1
q14=a1>b1
q15=a1<b1
q33=a3==b3
q32=a3>=b3
q33=a3<=b3
q34=a3>b3
q35=a3<b3

% declareGet(code,q01,'q01')
% declareGet(code,q02,'q02')
% declareGet(code,q03,'q03')
% declareGet(code,q04,'q04')
% declareGet(code,q05,'q05')
% declareGet(code,q11,'q11')
% declareGet(code,q12,'q12')
% declareGet(code,q13,'q13')
% declareGet(code,q14,'q14')
% declareGet(code,q15,'q15')
% declareGet(code,q31,'q31')
% declareGet(code,q32,'q32')
% declareGet(code,q33,'q33')
% declareGet(code,q34,'q34')
% declareGet(code,q35,'q35')

disp('tprod')
w0=tprod(a1,[-1],a3,[-1,1,2])
w1=tprod(w0,[1,2],a3,[3,1,2])
w11=tprod(w0,[1,2],a3,[3,1,2],'associate')

declareGet(code,w0,'w0')
declareGet(code,w1,'w1')
declareGet(code,w11,'w11')


disp('times, mtimes, rdivide, ldivide, rmdivide, rldivide')
z0=a0*b0
z01=a0*a0
z1=a0*a1
z2=a2*a1
z3=a2'*z2
z4=tprod(z3,[1])
z4=tprod(z3,[1],'associate')
z5=a3.*b3
z6=times(a3,b3,0)
z7=a3./b3
z8=a3.\b3
z9=(a2'*a2)\b1
z10=b1/(a2'*a2)

declareGet(code,z0,'z0')
declareGet(code,z01,'z01')
declareGet(code,z1,'z1')
declareGet(code,z2,'z2')
declareGet(code,z3,'z3')
declareGet(code,z4,'z4')
declareGet(code,z5,'z5')
declareGet(code,z6,'z6')
declareGet(code,z7,'z7')
declareGet(code,z8,'z8')
%declareGet(code,z9,'z9')
%declareGet(code,z10,'z10')

disp('min, max, all, any, abs')
r00=min(a0,[])
r01=min(a1,1)
r02=min(a3,[2,3])
r10=max(a0,[])
r11=max(a1,1)
r12=max(a3,[2,3])
r20=all(a0,[])
r21=all(a1,1)
r22=all(a3,[2,3])
r30=any(a0,[])
r31=any(a1,1)
r32=any(a3,[2,3]);
r4=abs(a0)
r5=abs(a2)

declareGet(code,r00,'r00')
declareGet(code,r01,'r01')
declareGet(code,r02,'r02')
declareGet(code,r10,'r10')
declareGet(code,r11,'r11')
declareGet(code,r12,'r12')
declareGet(code,r20,'r20')
declareGet(code,r21,'r21')
declareGet(code,r22,'r22')
declareGet(code,r30,'r30')
declareGet(code,r31,'r31')
declareGet(code,r32,'r32')
declareGet(code,r4,'r4')
declareGet(code,r5,'r5')

disp('diag, transpose')

s0=diag(a1)
s01=diag(a1,1)
s1=diag(a2*a2')
s10=diag(a2*a2',1)
s2=a0'
s3=a2'
s4=(a2')'
s5=tprod(a1,[1],Teye([3,3]),[1,2])    % essencially diag(vector->matrix)
s6=tprod(b2,[1,2],Teye([7,7]),[1,2])  % just non zero diagonal
s7=tprod(s6,[1,-1])
s8=tprod(s6,[1,-1],'associate')

declareGet(code,s0,'s0')
declareGet(code,s01,'s01')
declareGet(code,s1,'s1')
declareGet(code,s10,'s10')
declareGet(code,s2,'s2')
declareGet(code,s3,'s3')
declareGet(code,s4,'s4')
declareGet(code,s5,'s5')
declareGet(code,s6,'s6')
declareGet(code,s7,'s7')
declareGet(code,s8,'s8')

disp('chol, lu, ldl, pptrs, det')
v0=chol(b2)
v01=chol(b2,'typical.subscripts','typical.subscripts')
v1=lu(b2)
v2=ldl(b2)
v3=pptrs(a2'*a2,a1)
v4=det(ldl(a2'*a2))
v5=det(ldl(b2))

declareGet(code,v0,'v0')
declareGet(code,v01,'v01')
declareGet(code,v1,'v1')
declareGet(code,v2,'v2')
%declareGet(code,v3,'v3')
declareGet(code,v4,'v4')
%declareGet(code,v5,'v5')

disp('compose')
u0=exp(a0)
u1=log(a1)
u2=1./a2
u3=cube(a0)
u4=sqr(a1)
u5=sqrt(a2)
u6=cos(a3)
u7=sin(a3)
u8=tan(a3)
u9=atan(a3)
u10=normpdf(a3)

declareGet(code,u0,'u0')
declareGet(code,u1,'u1')
declareGet(code,u2,'u2')
declareGet(code,u3,'u3')
declareGet(code,u4,'u4')
declareGet(code,u5,'u5')
declareGet(code,u6,'u6')
declareGet(code,u7,'u7')
declareGet(code,u8,'u8')
declareGet(code,u9,'u9')
declareGet(code,u10,'u10')

disp('cat')
x1=[a1;b1]
x10=[a1;b1;7]
x2=[a2,b2]
%x20=[a2,b2,5]; % syntax error

declareGet(code,x1,'x1')
declareGet(code,x10,'x10')
declareGet(code,x2,'x2')
%declareGet(code,x20,'x20')

disp('gradient')
e0=a0
g0=gradient(e0,a0)
e01=5*a0+3
g01=gradient(e01,a0)
e02=7*a0*a0
g02=gradient(e02,a0)
e03=7*a0*a0+5*a0+3
g03=gradient(e03,a0)

e1=a1
g1=gradient(e1,a1)
e11=5*a1+3
g11=gradient(e11,a1)
e12=7*a1.*a1
g12=gradient(e12,a1)
e13=7*a1.*a1-5*a1+3
g13=gradient(e13,a1)

e2=c1*b2*c1+Tones(7)*c1
g2=gradient(e2,c1)
h2=gradient(g2,c1)
e21=tprod(c1,[-1],b2,[-1,-2],c1,[-2])+Tones(7)*c1
g21=gradient(e21,c1)
h21=gradient(g21,c1)

declareGet(code,e0,'e0')
declareGet(code,g0,'g0')
declareGet(code,e01,'e01')
declareGet(code,g01,'g01')
declareGet(code,e02,'e02')
declareGet(code,g02,'g02')
declareGet(code,e03,'e03')
declareGet(code,g03,'g03')

declareGet(code,e1,'e1')
declareGet(code,g1,'g1')
declareGet(code,e11,'e11')
declareGet(code,g11,'g11')
declareGet(code,e12,'e12')
declareGet(code,g12,'g12')
declareGet(code,e13,'e13')
declareGet(code,g13,'g13')

declareGet(code,e2,'e2')
declareGet(code,g2,'g2')
declareGet(code,h2,'h2')
declareGet(code,e21,'e21')
declareGet(code,g21,'g21')
declareGet(code,h21,'h21')


profile off;
profile viewer;