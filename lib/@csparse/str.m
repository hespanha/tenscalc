function outStr=str(obj,verboseLevel)
% outStr=str(obj,verboseLevel)
%   Produces a string representation of a csparse object.
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


    if nargin<2
        verboseLevel=0;
    end
    outStr='';
    if isempty(obj.sets) && isempty(obj.gets) && isempty(obj.copies) && isempty(obj.vectorizedOperations)
        outStr=[outStr,sprintf('  empty csparse object\n')];
        return;
    end
    if ~isempty(obj.vectorizedOperations)
        outStr=[outStr,sprintf('  elementary vectorized expressions (*=atomic):\n')];
        for i=1:height(obj.vectorizedOperations)
            type=getOne(obj.vectorizedOperations,'type',i);
            atomic=getOne(obj.vectorizedOperations,'atomic',i);
            description=getOne(obj.vectorizedOperations,'description',i);
            name=getOne(obj.vectorizedOperations,'name',i);
            osize=getOne(obj.vectorizedOperations,'osize',i);
            operands=getOne(obj.vectorizedOperations,'operands',i);
            subscripts=getOne(obj.vectorizedOperations,'subscripts',i);
            instructions=getOne(obj.vectorizedOperations,'instructions',i);
            if atomic
                outStr=[outStr,'*'];
            else
                outStr=[outStr,' '];
            end
            tc=find(obj.TCindex2CSvectorized==i);
            outStr=[outStr,sprintf('    %3d (TC:%-10s): %-20s %-12s = %20s(',i,index2str(tc),name,...
                             ['[',index2str(osize),']'],type)];
            ops=cellstr(getMulti(obj.vectorizedOperations,'name',operands));
            ops=sprintf(',%s',ops{:});
            outStr=[outStr,sprintf('%s)',ops(2:end))];

            if ~any(isnan(subscripts))
                if isempty(subscripts)
                    outStr=[outStr,sprintf('  nnz=%4d, fillin=%6.2f%%, size(subs)=[%s], size(instr)=[%s]',...
                                           size(subscripts,2),...
                                           100*size(subscripts,2)/prod(osize),...
                                           index2str(size(subscripts)),...
                                           index2str(size(instructions)))];
                else
                    if isempty(obj.memoryLocations)
                        outStr=[outStr,sprintf('  nnz=%4d, fillin=%6.2f%%, size(subs)=[%s], size(instr)=[%s], instr=[%d...]',...
                                               size(subscripts,2),...
                                               100*size(subscripts,2)/prod(osize),...
                                               index2str(size(subscripts)),...
                                               index2str(size(instructions)),...
                                               instructions(1))];
                    else
                        outStr=[outStr,sprintf('  nnz=%4d, fillin=%6.2f%%, size(subs)=[%s], size(instr)=[%s], instr=[%d/%d...]',...
                                               size(subscripts,2),...
                                               100*size(subscripts,2)/prod(osize),...
                                               index2str(size(subscripts)),...
                                               index2str(size(instructions)),...
                                               instructions(1),obj.memoryLocations(instructions(1)))];
                    end
                end
            elseif ~any(isnan(instructions)) && ~isempty(instructions)
                if isempty(obj.memoryLocations)
                    outStr=[outStr,sprintf('  size(instr)=[%s], instr=[%d...]',...
                                           index2str(size(instructions)),...
                                           instructions(1))];
                else
                    outStr=[outStr,sprintf('  size(instr)=[%s], instr=[%d/%d...]',...
                                           index2str(size(instructions)),...
                                           instructions(1),...
                                           obj.memoryLocations(instructions(1)))];
                end
            end
            if ~any(isnan(operands))
                outStr=[outStr,sprintf(', size(ops)=[%s]',index2str(size(operands)))];
            end
            outStr=[outStr,sprintf('\t%s\n',description)];
        end
    end
    if ~isempty(obj.sets)
        outStr=[outStr,sprintf('  sets:\n')];
        for i=1:length(obj.sets)
            outStr=[outStr,sprintf('    %20s(%10s(%2d)=...)          ',obj.sets(i).functionName,...
                             getOne(obj.vectorizedOperations,'name',obj.sets(i).destination),...
                             obj.sets(i).destination)];
            if ~isempty(obj.sets(i).childrenGroups)
                outStr=[outStr,sprintf(', children groups = %-20s',index2str(obj.sets(i).childrenGroups))];
            end
            outStr=[outStr,sprintf('\n')];
        end
    end
    if ~isempty(obj.gets)
        outStr=[outStr,sprintf('  gets:\n')];
        for i=1:length(obj.gets)
            outStr=[outStr,sprintf('    %20s(...=%10s(%2d))          ',obj.gets(i).functionName,...
                             getOne(obj.vectorizedOperations,'name',obj.gets(i).source),...
                             obj.gets(i).source)];
            if ~isempty(obj.gets(i).parentGroups)
                outStr=[outStr,sprintf(', parent groups = %-20s',index2str(obj.gets(i).parentGroups))];
            end
            outStr=[outStr,sprintf('\n')];
        end
    end
    if ~isempty(obj.copies)
        outStr=[outStr,sprintf('  copies:\n')];
        for i=1:length(obj.copies)
            outStr=[outStr,sprintf('    %20s(',obj.copies(i).functionName)];
            for j=1:length(obj.copies(i).destination)
                outStr=[outStr,sprintf('    %10s(%2d)=%10s(%2d),',...
                                 getOne(obj.vectorizedOperations,'name',obj.copies(i).destination(j)),...
                                 obj.copies(i).destination(j),...
                                 getOne(obj.vectorizedOperations,'name',obj.copies(i).source(j)),...
                                 obj.copies(i).source(j))];
            end
            outStr=[outStr(1:end-1),')'];
            if ~isempty(obj.copies(i).childrenGroups)
                outStr=[outStr,sprintf(', children groups = %-20s',index2str(obj.copies(i).childrenGroups))];
            end
            if ~isempty(obj.copies(i).parentGroups)
                outStr=[outStr,sprintf(', parent groups = %-20s',index2str(obj.copies(i).parentGroups))];
            end
            outStr=[outStr,sprintf('\n')];
        end
    end

    outStr=[outStr,sprintf('  instructions: %d instructions, %dx%d dependency graph with %d non-zero entries)\n',instructionsTableHeight(),size(obj.dependencyGraph),nnz(obj.dependencyGraph))];
    if verboseLevel>0
        for i=1:instructionsTableHeight()
            [type,parameters,operands]=getInstruction(int64(i));

            if strcmp(obj.typeInstructions,'matlab') && ~isempty(parameters)
                parameters=getArrayFromByteStream(uint8(parameters));
                if ~isnumeric(parameters)
                    parameters=nan;
                end
            end
            outStr=[outStr,sprintf('%8d/%8d: %d par=%-8s oprs=%-40s ',i,obj.memoryLocations(i),...
                             type,...
                             mymat2str(parameters),...
                             index2str(operands))];
            group=obj.instructionsGroup(i);
            if ~isnan(group)
                outStr=[outStr,sprintf('group=%d ',group)];
            end
            outStr=[outStr,sprintf('\n')];
        end
    end
end
