function computeMatlabInstructions(obj,ks)
% computeInstructions(obj,ks)
%   1) Computes the subscripts of the nonzero elements of the
%      vectorizedOperations indiced by ks (or all of them is k is
%      omitted)

%   2) Computes intructions to perform the computations needed for the
%      non-zero elements.
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

%obj.vectorizedOperations

    if nargin<2
        ks=1:height(obj.vectorizedOperations);
    end

    fast=false; % search for duplicate instructions

    %fprintf('  computeMatlabInstructions... ');
    %t0=clock;
    for thisExp=ks(:)'
        name=getOne(obj.vectorizedOperations,'name',thisExp);
        type=getOne(obj.vectorizedOperations,'type',thisExp);
        osize=getOne(obj.vectorizedOperations,'osize',thisExp);
        operands=getOne(obj.vectorizedOperations,'operands',thisExp);
        instructions=getOne(obj.vectorizedOperations,'instructions',thisExp);
        if ~any(isnan(instructions))
            continue;
        end

        if ~isempty(operands)
            computeMatlabInstructions(obj,operands);
        end

        parameters=getOne(obj.vectorizedOperations,'parameters',thisExp);
        opinstr=nan(length(operands),1);
        for i=1:length(operands)
            opinstr(i)=getOne(obj.vectorizedOperations,'instructions',operands(i));
        end

        switch type
          case {'variable'}
            if isempty(operands)
                % new variable
                instructions=newInstructions(obj,obj.Itypes.I_set,{[]},{[]},thisExp,fast);
            else
                % alias
                instructions=getOne(obj.vectorizedOperations,'instructions',operands);
            end

          case {'constant'}
            instructions=newInstructions(obj,obj.Itypes.I_load,...
                                         {double(getByteStreamFromArray({osize,parameters}))},...
                                         {[]},...
                                         thisExp,fast);

          otherwise
            Itype=getfield(obj.Itypes,sprintf('I_M%s',type));
            instructions=newInstructions(obj,Itype,...
                                         {double(getByteStreamFromArray({osize,parameters}))},...
                                         {opinstr},...
                                         thisExp,fast);
        end

        %fprintf('%3d: %s (instr=%d)\n',thisExp,name,instructions);

        set(obj.vectorizedOperations,'instructions',thisExp,instructions);
    end

    %fprintf('    done computeMatlabInstructions (%.3f sec)\n',etime(clock(),t0));

end
