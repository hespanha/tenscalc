function [subsY,instrY]=sparsity_subsref(obj,thisExp)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'subsref' and returns the positions in memory
%   of the nonzero elements
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

    verboseLevel=0;

    operands=getOne(obj.vectorizedOperations,'operands',thisExp);
    S=getOne(obj.vectorizedOperations,'parameters',thisExp);
    oosize=getOne(obj.vectorizedOperations,'osize',operands(1));
    subsX=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instrX=getOne(obj.vectorizedOperations,'instructions',operands(1));

    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('\n  sparsify_subsref(%3d): X  ndim=%d, nnz=%d\n',...
                thisExp,size(subsX,1),length(instrX));
    end

    switch S.type
      case '()',
        if length(S.subs)~=size(subsX,1)
            error(['mismatch between object length (%d) and ' ...
                   'indexing length (%d)'],size(subsX,1),length(S.subs));
        end
        %S.subs{:}
        ind_y=(1:length(S.subs{1}))';
        ind_subsref=S.subs{1}(:);
        for i=2:length(S.subs)
            % prune uninteresting rows in ind_y for next expansion
            k1=ismember(ind_subsref,subsX(1:i-1,:)','rows');
            ind_subsref=ind_subsref(k1,:);
            ind_y=ind_y(k1,:);

            % expand next dimension
            ind_subsref=[kron(ones(length(S.subs{i}),1),ind_subsref),...
                         kron(S.subs{i}(:),ones(size(ind_subsref,1),1))];
            ind_y   =[kron(ones(length(S.subs{i}),1),ind_y),...
                      kron((1:length(S.subs{i}))',ones(size(ind_y,1),1))];
        end
        % ind_subsref(1:min(20,end),:)
        % size(ind_subsref)
        % ind_y(1:min(20,end),:)
        % size(ind_y)
        % size(subsX)

        %tic
        [lia,k2]=ismember(ind_subsref,subsX','rows');
        k1=(1:size(ind_subsref,1))';
        k1=k1(lia);
        k2=k2(lia);
        % toc;
        % if any(any(ind_subsref(k1,:)'~=subsX(:,k2)))
        %     error('ismember mismatch');
        % end

        subsY=ind_y(k1,:)';
        instrY=instrX(k2);

      otherwise,
        error('subsref of type ''%s'' not implemented\n',S.type);
    end

    if verboseLevel>0
        fprintf('  sparsify_subsref(%3d): Y  ndim=%d, nnz=%d,  new instr=%4d (%d..%d) (%.2f sec)\n',...
        thisExp,size(subsY,1),length(instrY),...
            instructionsTableHeight()-n0,n0+1,instructionsTableHeight(),etime(clock,t0));
    end

end
