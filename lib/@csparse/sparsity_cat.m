function [subsY,instrY,instructions]=sparsity_cat(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'cat', allocates the required memory, and
%   determines the instructions needed to compute each of its nonzero
%   elements.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    osize=getOne(obj.vectorizedOperations,'osize',thisExp);
    operands=getOne(obj.vectorizedOperations,'operands',thisExp);
    dimension=getOne(obj.vectorizedOperations,'parameters',thisExp);

    %% Compute sparsity pattern
    subsY=zeros(length(osize),0,'uint64');
    instrY=zeros(0,1);

    j=0;
    for i=1:length(operands)
        subsX=getOne(obj.vectorizedOperations,'subscripts',operands(i));
        instrX=getOne(obj.vectorizedOperations,'instructions',operands(i));

        subsY=[subsY,[subsX(1:dimension-1,:);subsX(dimension,:)+j;subsX(dimension+1:end,:)]];

        osize1=getOne(obj.vectorizedOperations,'osize',operands(i));
        j=j+osize1(dimension);

        instrY=[instrY;instrX];
    end

    % if isequal(osize,[16,40])
    %     error('cat G_z')
    % end

    % if isequal(osize,[40,40])
    %     error('cat WW')
    % end
end
