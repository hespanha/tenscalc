function computeScalarInstructions(obj,ks,folder)
% computeInstructions(obj,ks)
%   1) Computes the subscripts of the nonzero elements of the
%      vectorizedOperations indiced by ks (or all of them is k is
%      omitted)

%   2) Computes intructions to perform the computations needed for the
%      non-zero elements.
%
% Any information about the instructions will be written in the
% given folder
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

%if nargin<2
%    ks=1:height(obj.vectorizedOperations);
%end

if nargin<3
    folder='.';
end

if length(ks)>=100
    fprintf('  computeScalarInstructions (%d expressions)...\n',length(ks));
end
t0=clock;
for jj=1:length(ks)
    thisExp=ks(jj);
    if mod(jj,100)==0
        fprintf('   computeScalarInstructions: currently at expr %d/%d...\n',jj,length(ks));
    end
    subscripts=getOne(obj.vectorizedOperations,'subscripts',thisExp);
    if ~any(isnan(subscripts))
         continue;
    end
    type=getOne(obj.vectorizedOperations,'type',thisExp);
    atomic=getOne(obj.vectorizedOperations,'atomic',thisExp);
    osize=getOne(obj.vectorizedOperations,'osize',thisExp);
    operands=getOne(obj.vectorizedOperations,'operands',thisExp);
    if ~isempty(operands)
        computeScalarInstructions(obj,operands);
    end
    %% check for atomic children
    atomicChildren=false;
    for i=1:length(operands)
        atomicChildren=getOne(obj.vectorizedOperations,'atomic',operands(i));
        if atomicChildren
            break;
        end
    end
    if atomic
        if atomicChildren
            switch type
              case 'variable'
                if isempty(operands)
                    % variable
                    error('computeScalarInstructions: type ''%s'' not allowed as atomic operator & children\n',type);
                else
                    % alias
                    subscripts=getOne(obj.vectorizedOperations,'subscripts',operands(1));
                    instructions=getOne(obj.vectorizedOperations,'instructions',operands(1));
                end
              otherwise
                error('computeScalarInstructions: type ''%s'' not allowed as atomic operator & children\n',type);
            end
        else % atomicChildren
            switch type
              case 'variable'
                if isempty(operands)
                    % variable
                    error('computeScalarInstructions: type ''%s'' not allowed as atomic operator\n',type);
                else
                    % alias
                    subscripts=getOne(obj.vectorizedOperations,'subscripts',operands(1));
                    instructions=getOne(obj.vectorizedOperations,'instructions',operands(1));
                end
              case {'lu','lu_sym'}
                subsX=getOne(obj.vectorizedOperations,'subscripts',operands(1));
                instrX=getOne(obj.vectorizedOperations,'instructions',operands(1));
                
                if length(osize)~=2 || osize(1)~=osize(2)
                    error('(atomic) lu operation only valid for square matrices (size=[%s])',...
                          index2str(osize));
                end
                
                %% compute compressed-column form
                [ji,k]=sortrows(subsX(end:-1:1,:)');
                % http://www.mathworks.com/help/matlab/apiref/mxsetir.html
                ir=ji(:,2)'-1; 
                % http://www.mathworks.com/help/matlab/apiref/mxsetjc.html
                jc=zeros(1,osize(2)+1);
                this=0;
                for l=1:osize(2)
                    this=this+sum(ji(:,1)==l);
                    jc(l+1)=this;
                end
                instrX=instrX(k); % order of jc;

                obj.atomicVariables(end+1).Ap=jc;
                obj.atomicVariables(end).Ai=ir;
                obj.atomicVariables(end).n=osize(2);

                %% add instruction
                subscripts=[];
                if strcmp(type,'lu')
                    itype=obj.Itypes.I_luS2A;
                else
                    itype=obj.Itypes.I_luS2Asym;
                end
                instructions=newInstruction(obj,itype,...
                                            double([length(obj.atomicVariables),osize(1)]),...
                                            instrX',...
                                            thisExp,obj.fastRedundancyCheck);
                obj.statistics.lu{end+1,1}=struct(...
                    'sizeA',osize,...
                    'nnzA',length(instrX),...
                    'uniqueNnzA',length(unique(instrX)));
              otherwise
                disp(obj)
                error('computeScalarInstructions: type ''%s'' not yet implemented as atomic operator & not atomic children\n',type);
            end % switch
        end % atomicChildren
    else % atomic
        if atomicChildren
            switch type
              case 'mldivide'
                subsLU=getOne(obj.vectorizedOperations,'subscripts',operands(1));
                instrLU=getOne(obj.vectorizedOperations,'instructions',operands(1));
                subsb=getOne(obj.vectorizedOperations,'subscripts',operands(2));
                instrb=getOne(obj.vectorizedOperations,'instructions',operands(2));
                if size(subsb,1)~=1
                    error('atomic mldivide(LU,b) only implemented for 1-d vector b, not size(subsb)=[%s]',index2str(size(subsb)));
                end
                if length(instrb)~=getOne(obj.vectorizedOperations,'osize',operands(2));
                    error('atomic mldivide(LU,b) only implemented for full b vector');
                end
                [subsb,k]=sort(subsb); % no need for sort rows since b has a single dimension
                instrb=instrb(k);

                subscripts=subsb;
                instructions=nan(length(instrb),1);
                % computes mldivide
                instructions(1)=newInstruction(obj,obj.Itypes.I_mldivideA2F1,...
                                               [],...
                                               [instrLU;instrb],...
                                               thisExp,obj.fastRedundancyCheck);
                % copy atomic mldivide result to scratchbook
                for i=2:length(instrb)
                    instructions(i)=newInstruction(obj,obj.Itypes.I_mldivideA2Fn,...
                                                     [i],...
                                                     [instrLU;instructions(1)],...
                                                     thisExp,obj.fastRedundancyCheck);
                end

                
              otherwise
                error('computeScalarInstructions: type ''%s'' not yet implemented as nonatomic operator & atomic children\n',type);
            end
        else % atomicChildren
            switch type
              case 'variable'
                if isempty(operands)
                    % variable
                    subscripts=memory2subscript(osize,1:prod(osize));
                    instructions=newInstructions(obj,obj.Itypes.I_set,...
                                                 cell(size(subscripts,2),1),{[]},...
                                                 thisExp,obj.fastRedundancyCheck);
                else
                    % alias
                    subscripts=getOne(obj.vectorizedOperations,'subscripts',operands(1));
                    instructions=getOne(obj.vectorizedOperations,'instructions',operands(1));
                end
              case 'constant'
                value=getOne(obj.vectorizedOperations,'parameters',thisExp);
                if length(osize)==2
                    % more efficient for 2 dimensions (especially if sparse)
                    [i,j,value]=find(value);
                    subscripts=[i(:),j(:)]';
                    value=value(:);
                else
                    value=value(:);   % flatten to 1st dimension
                    subscripts=memory2subscript(osize,1:prod(osize));
                    k=find(value==0); % remove zeros
                    if ~isempty(k)
                        subscripts(:,k)=[];
                        value(k)=[];
                    end
                end
                if ~isempty(value)
                    instructions=newInstructions(obj,obj.Itypes.I_load,...
                                                 num2cell(value),{[]},...
                                                 thisExp,obj.fastRedundancyCheck);
                else
                    instructions=[];
                    disp(obj.vectorizedOperations,thisExp)
                    fprintf('\n WARNING: computeScalarInstructions for empty constant (expression %d)\n',thisExp);
                end
                
              case 'full'
                subs1=getOne(obj.vectorizedOperations,'subscripts',operands(1));
                inst1=getOne(obj.vectorizedOperations,'instructions',operands(1));
                subscripts=memory2subscript(osize,1:prod(osize));
                [lia,locb]=ismember(subscripts',subs1','rows');
                instructions=nan(size(subscripts,2),1);
                instructions(lia)=inst1(locb(lia));
                if ~all(lia)
                    instructions(~lia)=newInstructions(obj,obj.Itypes.I_load,...
                                                       num2cell(zeros(sum(~lia),1)),{[]},...
                                                       thisExp,obj.fastRedundancyCheck);
                end
                    
              case 'zeros'
                subscripts=zeros(length(osize),0);
                instructions=[];
                
              case 'eye'
                osize1=osize(1:length(osize)/2);
                sub1=memory2subscript(osize1,1:prod(osize1));
                subscripts=[sub1;sub1];
                instructions=repmat(newInstructions(obj,obj.Itypes.I_load,...
                                                    {1},{[]},thisExp),prod(osize1),1);
                
              case 'ones'
                subscripts=memory2subscript(osize,1:prod(osize));
                instructions=repmat(newInstructions(obj,obj.Itypes.I_load,...
                                                    {1},{[]},thisExp),prod(osize),1);
                
              case 'reshape'
                subscripts=getOne(obj.vectorizedOperations,'subscripts',operands(1));
                oosize=getOne(obj.vectorizedOperations,'osize',operands(1));
                while length(oosize)<2
                    oosize(end+1)=1;
                    subscripts=[subscripts;ones(1,size(subscripts,2))];
                end
                subscripts=num2cell(subscripts,2);
                ind=sub2ind(oosize,subscripts{:});
                subscripts=memory2subscript(osize,ind);
                
                % no change in the instructions
                instructions=getOne(obj.vectorizedOperations,'instructions',operands(1));
                
              case 'vec2tensor'
                subscripts=getOne(obj.vectorizedOperations,'subscripts',operands(1));
                pars=getOne(obj.vectorizedOperations,'parameters',thisExp);
                dim=pars{1};
                sz=pars{2};
                subs=pars{3};
                subscripts=[subscripts(1:dim-1,:);
                            subs(:,subscripts(dim,:));
                            subscripts(dim+1:end,:)];
                % no change in the instructions
                instructions=getOne(obj.vectorizedOperations,'instructions',operands(1));                

                % keep subscripts in the "natural" order (sorted by row and then col)
                [subscripts,k]=sortrows(subscripts',size(subscripts,1):-1:1);
                subscripts=subscripts';
                instructions=instructions(k);
                
              case 'subsref'
                [subscripts,instructions]=sparsity_subsref(obj,thisExp);
                
              case 'cat'
                [subscripts,instructions]=sparsity_cat(obj,thisExp);
                
              case 'diag'
                [subscripts,instructions]=sparsity_diag(obj,thisExp);
                
              case 'min'
                [subscripts,instructions]=sparsity_min(obj,thisExp);
                
              case 'max'
                [subscripts,instructions]=sparsity_max(obj,thisExp);
                
              case 'min2'
                [subscripts,instructions]=sparsity_min2(obj,thisExp);
                
              case 'max2'
                [subscripts,instructions]=sparsity_max2(obj,thisExp);
                
              case 'plus'
                [subscripts,instructions]=sparsity_plus(obj,thisExp);
            
              case 'norm2'
                [subscripts,instructions]=sparsity_norm2(obj,thisExp);
            
              case 'norm1'
                [subscripts,instructions]=sparsity_norm1(obj,thisExp);
            
              case 'norminf'
                [subscripts,instructions]=sparsity_norminf(obj,thisExp);
            
              case 'clp'
                [subscripts,instructions]=sparsity_clp(obj,thisExp);
            
              case 'rdivide'
                [subscripts,instructions]=sparsity_rdivide(obj,thisExp);
            
              case 'tprod'
                [subscripts,instructions]=sparsity_tprod(obj,thisExp);
            
              case {'transpose','ctranspose'}
                [subscripts,instructions]=sparsity_transpose(obj,thisExp);
            
              case 'lu'
                files=getOne(obj.vectorizedOperations,'parameters',thisExp);
                files{3}=fullfile(folder,files{3});
                files{4}=fullfile(folder,files{4});
                [subscripts,instructions,p,q]=sparsity_lu(obj,thisExp,files{3},files{4});
                set(obj.vectorizedOperations,'parameters',thisExp,{p,q,files{3},files{4}});
            
              case 'ldl'
                files=getOne(obj.vectorizedOperations,'parameters',thisExp);
                files{3}=fullfile(folder,files{3});
                files{4}=fullfile(folder,files{4});
                [subscripts,instructions,p]=sparsity_ldl(obj,thisExp,files{3},files{4});
                set(obj.vectorizedOperations,'parameters',thisExp,{p,files{3},files{4}});
                
                % case 'chol'
                %   files=getOne(obj.vectorizedOperations,'parameters',thisExp);
                %   files{3}=fullfile(folder,files{3});
                %   files{4}=fullfile(folder,files{4});
                %   [subscripts,instructions,p]=sparsity_chol(obj,thisExp,files{3},files{4});
                %   set(obj.vectorizedOperations,'parameters',thisExp,{p,files{3},files{4}});
                
              case 'mldivide_l1'
                [subscripts,instructions]=sparsity_mldivide_l1(obj,thisExp);
                
              case 'mldivide_u1'
                [subscripts,instructions]=sparsity_mldivide_u1(obj,thisExp);
                
              case 'mldivide_u'
                [subscripts,instructions]=sparsity_mldivide_u(obj,thisExp);
                
              case 'mldivide_d'
                [subscripts,instructions]=sparsity_mldivide_d(obj,thisExp);
                
              case 'lu_d'
                [subscripts,instructions]=sparsity_lu_d(obj,thisExp);
                
              case 'lu_l'
                [subscripts,instructions]=sparsity_lu_l(obj,thisExp);
                
              case 'lu_u'
                [subscripts,instructions]=sparsity_lu_u(obj,thisExp);
                
              case 'ldl_d'
                [subscripts,instructions]=sparsity_ldl_d(obj,thisExp);
                
              case 'ldl_l'
                [subscripts,instructions]=sparsity_ldl_l(obj,thisExp);
                
              case 'logdet_ldl'
                [subscripts,instructions]=sparsity_logdet_ldl(obj,thisExp);
                
              case 'logdet_lu'
                [subscripts,instructions]=sparsity_logdet_lu(obj,thisExp);
                
              case 'det_ldl'
                [subscripts,instructions]=sparsity_det_ldl(obj,thisExp);
                
              case 'det_lu'
                [subscripts,instructions]=sparsity_det_lu(obj,thisExp);
                
              case 'componentwise'
                funs=getOne(obj.vectorizedOperations,'parameters',thisExp);
                [subscripts,instructions]=sparsity_componentwise(obj,thisExp,funs);
                
              case 'compose'
                functions={'@(x__)log(x__)';'@(x__)1./x__';'@(x__)-1./x__.^2';
                           '@(x__)2./x__.^3';
                           '@(x__)exp(x__)';
                           '@(x__)sin(x__)';'@(x__)-sin(x__)';
                           '@(x__)cos(x__)';'@(x__)-cos(x__)';
                           '@(x__)round(x__)';'@(x__)ceil(x__)';'@(x__)floor(x__)';
                           '@(x__)abs(x__)';'@(x__)sign(x__)';
                           '@(x__)sqrt(x__)';'@(x__).5./sqrt(x__)';'@(x__)-.25./x__.^1.5';
                           '@(x__)x__.^2';'@(x__)2*x__';'@(x__)2*ones(size(x__))';
                           '@(x__)x__.^3';'@(x__)3*x__.^2';'@(x__)6*x__';
                           '@(x__)atan(x__)';'@(x__)1./(1+x__.^2)';'@(x__)-2*x__./(1+x__.^2).^2';
                           '@(x__)log(1+exp(x__))';'@(x__)1./(1+exp(-x__))';...
                                '@(x__)1./(2+exp(-x__)+exp(x__))';
                           '@(x__)relu(x__)';'@(x__)heaviside(x__)';'@(x__)zeros(size(x__))'};
                instruction={obj.Itypes.I_log;obj.Itypes.I_inv;obj.Itypes.I_minus_inv_sqr;
                             obj.Itypes.I_2_inv_cube;
                             obj.Itypes.I_exp;
                             obj.Itypes.I_sin;obj.Itypes.I_minus_sin;
                             obj.Itypes.I_cos;obj.Itypes.I_minus_cos;
                             obj.Itypes.I_round;obj.Itypes.I_ceil;obj.Itypes.I_floor;
                             obj.Itypes.I_abs;obj.Itypes.I_sign;
                             obj.Itypes.I_sqrt;obj.Itypes.I_Dsqrt;obj.Itypes.I_DDsqrt;
                             obj.Itypes.I_sqr;obj.Itypes.I_2times;obj.Itypes.I_2;
                             obj.Itypes.I_cube;obj.Itypes.I_3sqr;obj.Itypes.I_6times;
                             obj.Itypes.I_atan;obj.Itypes.I_Datan;obj.Itypes.I_DDatan;
                             obj.Itypes.I_srelu;obj.Itypes.I_dsrelu;obj.Itypes.I_ddsrelu;
                             obj.Itypes.I_relu;obj.Itypes.I_heaviside;obj.Itypes.I_zero};
                sparsityType={@sparsity_compose_full;@sparsity_compose_full;@sparsity_compose_full;
                              @sparsity_compose_full;
                              @sparsity_compose_full;
                              @sparsity_compose;@sparsity_compose;
                              @sparsity_compose_full;@sparsity_compose_full;
                              @sparsity_compose;@sparsity_compose;@sparsity_compose;
                              @sparsity_compose;@sparsity_compose_full;
                              @sparsity_compose;@sparsity_compose_full;@sparsity_compose_full;
                              @sparsity_compose;@sparsity_compose;@sparsity_compose_full;
                              @sparsity_compose;@sparsity_compose;@sparsity_compose;
                              @sparsity_compose;@sparsity_compose_full;@sparsity_compose;
                              @sparsity_compose;@sparsity_compose_full;@sparsity_compose_full;
                              @sparsity_compose;@sparsity_compose;@sparsity_compose};
                
                if length(functions)~=length(instruction) || length(functions)~=length(sparsityType)
                    error('computeScalarInstructions: internal inconsistency detected in compose()')
                end
                fun=getOne(obj.vectorizedOperations,'parameters',thisExp);
                k=find(strcmp(fun,functions));
                if isempty(k)
                    error('computeScalarInstructions: compose object not implemented for function ''%s''\n',fun);
                end
                if instruction{k}~=obj.Itypes.I_zero
                    [subscripts,instructions]=sparsityType{k}(obj,thisExp,instruction{k});
                else
                    osize=getOne(obj.vectorizedOperations,'osize',thisExp);
                    subscripts=zeros(length(osize),0);
                    instructions=zeros(0,1);
                end
              otherwise
                disp(obj)
                error('computeScalarInstructions: type ''%s'' not yet implemented\n',type);
            end % switch

        if ~isequal(size(subscripts,1),length(osize))
            obj
            error('addTCExpression: internal error object (%s) size [%s] does not match size of subscripts [%s]\n',type,index2str(osize),index2str(size(subscripts)));
        end
      end % atomicChildren
    end % atomic
        
    if ~isequal(size(subscripts,2),length(instructions)) && ... % tensor of non-trivial size
            ~(isempty(subscripts)&&length(instructions)==1)          % scalar
        % disp(obj)
        subscripts
        instructions
        error('addTCExpression: internal error object (%s) size of subscripts [%s] does not match size of instructions [%s]\n',type,index2str(size(subscripts)),index2str(size(instructions)));
    end

    if any(isnan(instructions)) 
        subscripts
        instructions
        error('computeScalarInstructions: nan instructions for type ''%s''\n',type);
    end

    if any(isnan(subscripts(:))) 
        subscripts
        instructions
        error('computeScalarInstructions: nan subscripts for type ''%s''\n',type);
    end

    set(obj.vectorizedOperations,'subscripts',thisExp,subscripts);
    set(obj.vectorizedOperations,'instructions',thisExp,instructions);
end

%fprintf('    done computeScalarInstructions (%.3f sec)\n',etime(clock(),t0));

end
