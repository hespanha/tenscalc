function computeMatlabInstructions(obj,ks)
% computeInstructions(obj,ks)
%   1) Computes the subscripts of the nonzero elements of the
%      vectorizedOperations indiced by ks (or all of them is k is
%      omitted)

%   2) Computes intructions to perform the computations needed for the
%      non-zero elements.
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
