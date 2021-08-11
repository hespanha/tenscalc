function [subsY,instrY]=sparsity_transpose(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'subsref' and returns the positions in memory
%   of the nonzero elements
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    verboseLevel=0;
    
    operands=getOne(obj.vectorizedOperations,'operands',thisExp);
    
    subsX=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instrX=getOne(obj.vectorizedOperations,'instructions',operands(1));
    
    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('\n  sparsify_transpose(%3d): X  ndim=%d, nnz=%d\n',...
                thisExp,size(subsX,1),length(instrX));
    end
    
    if size(subsX,1)~=2
        error('transpose only valid for 2D matrices (%d dimensions instead)',size(subsX,1));
    end
    subsY=subsX([2,1],:);
    instrY=instrX;
    
    if verboseLevel>0
        fprintf('  sparsify_transpose(%3d): Y  ndim=%d, nnz=%d,  new instr=%4d (%d..%d) (%.2f sec)\n',...
                thisExp,size(subsY,1),length(instrY),...
                instructionsTableHeight()-n0,n0+1,instructionsTableHeight(),etime(clock,t0));
    end
end
