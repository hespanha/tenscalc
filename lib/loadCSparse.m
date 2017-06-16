function [osize,subscripts,values,matrix]=loadCSparse(filename_subscripts,filename_values,ignoremagic)
% [osize,subscripts,values,matrix]=loadCSparse(filename_subscripts,filename_values)
%
% Load a sparse tensor saved by @csparse generated code.
%
% Inputs:
%   
%   filename_subscripts - filename with the subscripts of the
%                         nonzero entries
%   filename_values     - filename with the values of the
%                         nonzero entries
%
% Outputs:
%
%   osize      - dimension of the sparse tensor
%   subscripts - int64-valued array with the subscripts of the 
%                nonzero entries, with one row per dimension and
%                one column per nonzero entry
%   values     - double-valued array with the values of the nonzero entries
%   matrix (optional) - sparse matrix (only possible if the tensor
%                       is two dimensional (matrix) 
%
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

if nargin<3
    ignoremagic=false;
end

osize=nan;
subscripts=nan;
values=nan;
matrix=nan;

fd1=fopen(filename_subscripts,'r');
if fd1<0
    fprintf('    loadCsparse: cannot open subscripts file "%s" (fd=%d)\n',filename_subscripts,fd1)
    return
end

fd2=fopen(filename_values,'r');
if fd2<0
    fprintf('    loadCsparse: cannot open values file "%s" (fd=%d)\n',filename_values,fd2)
    fclose(fd1);
    return
end

magic1=fread(fd1,1,'int64=>int64');
magic2=fread(fd2,1,'int64=>int64');

if ~ignoremagic && magic1~=magic2
    fprintf('    loadCsparse: mismatch between magic numbers in subscripts "%s" and values "%s" files (%d ~= %d)\n',filename_subscripts,filename_values,magic1,magic2)
    fclose(fd1);
    fclose(fd2);
    return
end

nDim=fread(fd1,1,'int32=>int32');
osize=fread(fd1,nDim,'int64=>int64');
nnz=fread(fd1,1,'int64=>int64');
subscripts=fread(fd1,[nDim,nnz],'int64=>int64');
values=fread(fd2,[nnz,1],'double=>double');

fclose(fd1);
fclose(fd2);

if size(subscripts,2)~=length(values)
    fprintf('    loadCsparse: mismatch between subscripts file "%s" and values file "%s"\n',...
        filename_subscripts,filename_values);
    osize=nan;
    return
end

if nargout>=3
    while length(osize)<2
        osize(end+1)=1;
        subscripts=[subscripts;ones(1,size(subscripts,2),'int64')];
    end
    if length(osize)==2
        matrix=sparse(double(subscripts(1,:)),double(subscripts(2,:)),values,double(osize(1)),double(osize(2)));
    else
        osize
        error('can only return sparse 2-dimensional matrices');
    end
end


