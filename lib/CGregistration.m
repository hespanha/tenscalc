function CG=CGregistration()

% Copyright 2012-2017 Joao Hespanha

% This file is part of Tencalc.
%
% TensCalc is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% TensCalc is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with TensCalc.  If not, see <http://www.gnu.org/licenses/>.

%% Magics
CG.magic(1,1) = struct('type',16,'name','int16l','description','indices stored as 16bit integers with the least significant byte first (little-endian byte ordering)');
CG.magic(2,1) = struct('type',32,'name','int32l','description','indices stored as 32bit integers with the least significant byte first (little-endian byte ordering)');
CG.magic(3,1) = struct('type',64,'name','int64l','description','indices stored as 64bit integers with the least significant byte first (little-endian byte ordering)');
    
%% Constants
CG.constant(1,1)    =struct('type',10001,'name','double_full_tensor',...
                            'description','tensor with double-values entries [# dimension,size vector,double values]');
CG.constant(end+1,1)=struct('type',10002,'name','double_sparse_tensor',...
                            'description','tensor with double-values entries [# dimension,size vector,nonzero subscripts,double values]');

CG.constant(end+1,1)=struct('type',20001,'name','int64l_full_vector',...
                            'description','vector with int64l-values entries [int64l values]');
CG.constant(end+1,1)=struct('type',20002,'name','int64l_full_tensor',...
                            'description','vector with int64l-values entries [int64l values]');

CG.constant(end+1,1)=struct('type',30001,'name','string',...
                            'description','array of chars [char values]');

%% Vector Functions (1st integer-valued arrays, 2nd sparse tensor parameters)
CG.function(1    ,1)=struct('type',00000,'name','constant'    ,'description','constant function');

CG.function(end+1,1)=struct('type',10001,'name','Vzeros'       ,'description','Tzeros(size)');
CG.function(end+1,1)=struct('type',10002,'name','Vones'        ,'description','Tones(size)');
CG.function(end+1,1)=struct('type',10003,'name','Veye'         ,'description','Teye(size)');

CG.function(end+1,1)=struct('type',10005,'name','Vrepmat'      ,'description','repmat(size,a)');
CG.function(end+1,1)=struct('type',10006,'name','Vcat'         ,'description','cat(dim,a,b,...)');
CG.function(end+1,1)=struct('type',10007,'name','Vreshape'     ,'description','reshape(size,a)');
CG.function(end+1,1)=struct('type',10029,'name','Vsum'         ,'description','sum(dim,a,b,...)');
CG.function(end+1,1)=struct('type',10004,'name','Vsubsref'     ,'description','subsref(i1,i2,...,a)');
CG.function(end+1,1)=struct('type',10014,'name','Vmin'         ,'description','min(dim,a,b,...)');
CG.function(end+1,1)=struct('type',10015,'name','Vmax'         ,'description','max(dim,a,b,...)');
CG.function(end+1,1)=struct('type',10016,'name','Vall'         ,'description','all(dim,a,b,...)');
CG.function(end+1,1)=struct('type',10017,'name','Vany'         ,'description','any(dim,a,b,...)');

CG.function(end+1,1)=struct('type',10008,'name','Vfull'        ,'description','full(a)');
CG.function(end+1,1)=struct('type',10012,'name','Vnorm2'       ,'description','norm2(a)');
CG.function(end+1,1)=struct('type',10013,'name','Vnorminf'     ,'description','norminf(a)');
CG.function(end+1,1)=struct('type',10026,'name','Vctranspose'  ,'description','ctranspose(a)');
CG.function(end+1,1)=struct('type',10027,'name','Vtranspose'   ,'description','transpose(a)');
CG.function(end+1,1)=struct('type',10018,'name','Vabs'         ,'description','abs(a)');
CG.function(end+1,1)=struct('type',10021,'name','Vlu'          ,'description','lu(a)');
CG.function(end+1,1)=struct('type',10022,'name','Vchol'        ,'description','chol(a)');
CG.function(end+1,1)=struct('type',10030,'name','Vdiag'        ,'description','diag(a)');
CG.function(end+1,1)=struct('type',10023,'name','Vmldivide_l1' ,'description','mldivide_l1(a,b)');
CG.function(end+1,1)=struct('type',10024,'name','Vmldivide_u'  ,'description','mldivide_u(a,b)');
CG.function(end+1,1)=struct('type',10019,'name','Vclp'         ,'description','clp(a,b)');
CG.function(end+1,1)=struct('type',10028,'name','Vrdivide'     ,'description','rdivide(a,b)');

CG.function(end+1,1)=struct('type',10009,'name','Vplus'        ,'description','plus(signs,a,b,...)');
CG.function(end+1,1)=struct('type',10010,'name','Vtprod'       ,'description','tprod(size,sumsize,a,ia,b,ib,...)');

CG.function(end+1,1)=struct('type',10020,'name','Vcompose'     ,'description','compose(fun,a)');

CG.function(end+1,1)=struct('type',10011,'name','Vtprod_matlab','description','tprod_matlab(...)');
CG.function(end+1,1)=struct('type',10025,'name','Vmtimes'      ,'description','mtimes(...)');

%% Scalar Functions
CG.function(end+1,1)=struct('type',20001,'name','Ssum','description',...
                            'sum(parameters.*operands), parameters = vector with +1/-1 specifying addition/subtraction for each element of vector');

CG.function(end+1,1)=struct('type',20002,'name','Ssumprod','description',...
        'sum(prod(operands,1))=a1*a2*...*am+b1*b2*...*bm+..., operands = a1,a2,...,am,b1,b2,...,bm,c1,... In case operands are one long vector: parameters = [m; # sums]');
    
CG.function(end+1,1)=struct('type',20003,'name','Sinv','description',...
                            'inverse: 1/denominator, operands = denominator');
CG.function(end+1,1)=struct('type',20004,'name','Sminus_inv_sqr','description',...
                            'inverse: -1/denominator^2, operands = denominator');
CG.function(end+1,1)=struct('type',20005,'name','Sdiv','description',...
                            'division: numerator/denominator, operands = numerator, denominator');
CG.function(end+1,1)=struct('type',20006,'name','Sminus_dot','description',...
                            'symmetric of dot-product: -a1*b1-a2*b2-..., operands = a1,b1,a2,b2,....');
CG.function(end+1,1)=struct('type',20007,'name','Sminus_dot_div','description',...
                            'symmetric of dot-product & divide: (-a1*b1-a2*b2-...)/a0, operands = a0,a1,b1,a2,b2,....');
CG.function(end+1,1)=struct('type',20008,'name','Splus_minus_dot','description',...
                            'subtract dot-product: a0-a1*b1-a2*b2-..., operands = a0,a1,b1,a2,b2,....');
CG.function(end+1,1)=struct('type',20009,'name','Splus_sqr','description',...
                            'norm-2 square: a1^2+a2^2+..., operands = a1,a2,....');
CG.function(end+1,1)=struct('type',20010,'name','Splus_minus_dot_div',...
                            'description','subtract dot-product & divide: (a0-a1*b1-a2*b2-...)/b0, operands = a0,b0,a1,b1,a2,b2,....');
CG.function(end+1,1)=struct('type',20011,'name','Smin','description',...
                            'computes the min: min(a1,a2,...);, operands = a1,a2,....');
CG.function(end+1,1)=struct('type',20012,'name','Smin0','description',...
                            'computes the min: min(0,a1,a2,...);, operands = a1,a2,....');
CG.function(end+1,1)=struct('type',20013,'name','Smax','description',...
                            'computes the max: max(a1,a2,...);, operands = a1,a2,....');
CG.function(end+1,1)=struct('type',20014,'name','Smax0','description',...
                            'computes the max: max(0,a1,a2,...);, operands = a1,a2,....');
CG.function(end+1,1)=struct('type',20015,'name','Smax_abs','description',...
                            'infinity norm: max{ |a1|,|a2|,... }, operands = a1,a2,....');
CG.function(end+1,1)=struct('type',20016,'name','Sclp','description',...
                            'computes the largest c>=0 such that: ai + c bi >=0, operands = a1,b1,a2,b2,....an,bn');
CG.function(end+1,1)=struct('type',20017,'name','Sexp','description',...
                            'computes exp(a1), operand = a1');
CG.function(end+1,1)=struct('type',20018,'name','Slog','description',...
                            'computes log(a1), operand = a1');
CG.function(end+1,1)=struct('type',20019,'name','Scos','description',...
                            'computes cos(a1), operand = a1');
CG.function(end+1,1)=struct('type',20020,'name','Sminus_cos','description',...
                            'computes -cos(a1), operand = a1');
CG.function(end+1,1)=struct('type',20021,'name','Ssin','description',...
                            'computes sin(a1), operand = a1');
CG.function(end+1,1)=struct('type',20022,'name','Sminus_sin','description',...
                            'computes -sin(a1), operand = a1');
CG.function(end+1,1)=struct('type',20023,'name','Sabs','description',...
                            'computes abs(a1), operand = a1');
CG.function(end+1,1)=struct('type',20024,'name','Ssqrt','description',...
                            'computes sqrt(a1), operand = a1');
CG.function(end+1,1)=struct('type',20025,'name','SDsqrt','description',...
                            'computes .5/sqrt(a1), operand = a1');
CG.function(end+1,1)=struct('type',20026,'name','SDDsqrt','description',...
                            'computes -.25/sqrt(a1)^3, operand = a1');
CG.function(end+1,1)=struct('type',20027,'name','Ssqr','description',...
                            'computes a1^2, operand = a1');
CG.function(end+1,1)=struct('type',20028,'name','S2times','description',...
                            'computes 2*a1, operand = a1');
CG.function(end+1,1)=struct('type',20029,'name','S2','description',...
                            'computes 2, operand = a1 (ignored)');
CG.function(end+1,1)=struct('type',20030,'name','Satan','description',...
                            'computes atan(a1), operand = a1');
CG.function(end+1,1)=struct('type',20031,'name','SDatan','description',...
                            'computes 1/(1+a1^2), operand = a1');
CG.function(end+1,1)=struct('type',20032,'name','SDDatan','description',...
                            'computes -2*a1/(1+a1^2)^2, operand = a1');
    
%% tensor-scalar convertion functions
CG.function(end+1,1)=struct('type',30001,'name','scalars2ftensor','description',...
                            'compiles scalars into a full variable: tensor=f(subscripts,a1,a2,a3,...)');
CG.function(end+1,1)=struct('type',30002,'name','ftensor2scalars','description',...
                            'compiles scalars into a full variable: [a1,a2,a3,...]=f(tensor)');
    
%% I/O functions
CG.io(1,1)    =struct('type', 1,'name','set'    ,'description','set');
CG.io(end+1,1)=struct('type', 2,'name','get'    ,'description','get');
CG.io(end+1,1)=struct('type', 3,'name','copy'   ,'description','copy');

%% Symbols
CG.symbol(1,1)    =struct('type', 1,'name','constant','description','constant symbol');
CG.symbol(end+1,1)=struct('type', 2,'name','function','description','function-node symbol');
CG.symbol(end+1,1)=struct('type', 3,'name','variable','description','variable-node symbol');
CG.symbol(end+1,1)=struct('type', 4,'name','io',      'description','io-function symbol');

%% Save as json
%savejson('',CG,struct('FileName','CGregistration.json','SingleArray',1));
end