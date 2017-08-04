function [subsLDL,instrLDL,p]=sparsity_ldl(obj,thisExp,typical_subscripts,typical_values)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'ldl' and returns the positions in memory
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

    % profile on

    verboseLevel=2;  
    
    operands=getOne(obj.vectorizedOperations,'operands',thisExp);
    
    subsX=getOne(obj.vectorizedOperations,'subscripts',operands(1));
    instrX=getOne(obj.vectorizedOperations,'instructions',operands(1));
    % typeX=getMulti(obj.instructions,'type',instrX);
    % kSet=typeX==obj.Itypes.I_set;
    % kLoad=typeX==obj.Itypes.I_load';
    osize=getOne(obj.vectorizedOperations,'osize',thisExp);
    n=osize(1);
    
    if length(osize)~=2 || osize(1)~=osize(2)
        error('LDL decomposition only implemented for square matrices\n');
    end
    
    [uniqueInstrX,~,kuniqueInstrX]=unique(instrX); % instrX=uniqueInstrX(kuniqueInstrX);

    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsifyLDL(%3d): A  size=%-10s, nnz=%d, distinct nnz=%d\n',...
                thisExp,['[',index2str(osize),']'],length(instrX),length(uniqueInstrX));
    end

    
    % compute matrix to optimize sparsity in ldl factorization
    [sz,subscripts,values,A]=loadCSparse(typical_subscripts,typical_values);

    if any(isnan(sz)) || ~isequal(subscripts,subsX) 
        fprintf('    using random values, ');

        s = RandStream('mt19937ar','Seed',0);
        % create random values preserving structural symmetries
        values=1+.5*rand(s,length(uniqueInstrX),1);
        values=values(kuniqueInstrX);
        
        A=sparse(double(subsX(1,:)),double(subsX(2,:)),values,n,n);
    else
        fprintf('    using values from "%s", ',typical_values);
    end


    nNonSymA=symmetric(obj,subsX,instrX);
    if nNonSymA==0
        fprintf('structurally symmetric matrix, ');
        if ~isequal(A,A')
            norm(full(A-A'))
            error('Unexpected non symmetric matrix');
        end
    else
        fprintf('%d/%d nonzero entries of A and A'' DO NOT MATCH - LOWER ONES USED\n',...
                nNonSymA,length(instrX));
    end
    
    % determine column permutation for "optimal sparsity"
    
    if 0
        % sometimes does not guarantee existance of LDL with D
        % diagonal (instead of block diagonal)
        [lp,dp,p]=ldl((A+A')/2,.01,'vector');  % .01 gives more flexibility in minimizing fillin
    elseif 0
        % sometimes gets very bad filling
        [lp,up,p]=lu((A+A')/2,.01,'vector');  % .01 gives more
                                              % flexibility in
                                              % minimizing fillin;
    elseif 1
        % larger threshold (default 10) for dense (slower, but more accurate)
        % smvm(30): 27335 nnz A/15016 unique nnz L/102766 unique instructions
        p=symamd((A+A')/2,[1000,0]); 
        %p=symamd((A+A')/2,[1000,1]); % verbose
    elseif 0
        % faster than symamd, but seems to lead to a little more fillin
        % smvm(30): 27335 nnz A/15646 unique nnz L/104026 unique instructions
        % larger threshold (default 10) for dense (slower, but more accurate)
        p=amd((A+A')/2,struct('dense',1000,'aggressive',true));
    end
        
    if verboseLevel>1
        fig=gcf;
        f=figure(101);
        set(f,'Name',typical_subscripts);
        clf

        subplot(2,3,1)
        spy(A);
        if verboseLevel>2
            fprintf('\nA\n');
            disp(full(A(1:10,1:10)))
        end
        % plot(subsX(1,kSet),subsX(2,kSet),'r+')
        % plot(subsX(1,kLoad),subsX(2,kLoad),'g+')
        if nNonSymA==0
            title(sprintf('symmetric A %dx%d (original)',size(A)));
        else
            title(sprintf('asymetric A %dx%d (original)',size(A)));
        end
        if exist('lp','var')
            subplot(2,3,2)
            spy(abs(lp));
            title('L (matlab''s ld with perm)');
            if verboseLevel>2
                fprintf('L (matlab''s ld with perm\n')
                disp(full(lp(1:10,1:10)))
            end
        end

        subplot(2,3,3)
        plot(1:length(p),p,'rx');
        legend('row/col  perm');
        axis ij;axis image;
        
        subplot(2,3,4)
        spy(A(p,p))
        title('A (permuted)');
        if verboseLevel>2
            fprintf('A (permuted)\n');
            disp(full(A(p(1:10),p(1:10))))
        end
        drawnow
    end
 
    %% Compute instructions
    
    if any(instrX==0)
        error('sparsityLDL: algorithm expects no instrX=0\n');
    end
    
    % permute X
    instrXp=sparse(double(subsX(1,:)),...
                   double(subsX(2,:)),...
                   double(instrX),n,n);
    instrXp=instrXp(p,p);
    %full(instrXp)

    % initialize LDL
    instrLDL=sparse([],[],[],n,n);
    instrD=zeros(1,n);
    
    for col=1:n
        % v=nan(col-1,1);
        % for i=1:col-1 
        %     v(i)=LDL(col,i)*D(i);
        % end 
        instrV=zeros(1,col-1);
        is=find(instrLDL(col,1:col-1) & instrD(1:col-1));
        if ~isempty(is)
            parameters=repmat({[2,1]},1,length(is));
            operands=num2cell(full([instrLDL(col,is);instrD(is)]),1);
            instrV(is)=newInstructions(obj,obj.Itypes.I_sumprod,repmat({[2,1]},1,length(is)),...
                                       operands,...
                                       thisExp,obj.fastRedundancyCheck);
        end
        % D(col)=A(col,col)- LDL(col,1:col-1)*v(1:col-1);
        is=find(instrLDL(col,1:col-1) & instrV);
        if instrXp(col,col)
            if isempty(is)
                instrD(col)=instrXp(col,col);
            else
                operands=full([instrXp(col,col),reshape([instrLDL(col,is);instrV(is)],1,[])]);
                instrD(col)=newInstruction(obj,obj.Itypes.I_plus_minus_dot,[],...
                                           operands,...
                                           thisExp,obj.fastRedundancyCheck);
            end
        else
            if isempty(is)
                error('ldl needs pivoting');
            end
            operands=reshape(full([instrLDL(col,is);instrV(is)]),1,[]);
            instrD(col)=newInstruction(obj,obj.Itypes.I_minus_dot,[],...
                                       operands,...
                                       thisExp,obj.fastRedundancyCheck);
        end
        
        
        % for row=col+1:n
        %     LDL(row,col)=( A(row,col)- LDL(row,1:col-1)*v ) / D(col) ;
        %     LDL(col,row)=LDL(row,col); % not needed, but just to get L'
        % end

        rows=find( instrXp(col+1:n,col) | instrLDL(col+1:n,1:col-1)*instrV' )';
        for row=col+rows
            is=find(instrLDL(row,1:col-1) & instrV);
            if instrXp(row,col)
                if isempty(is) 
                    operands=full([instrXp(row,col),instrD(col)]);
                    instrLDL(row,col)=newInstruction(obj,obj.Itypes.I_div,[],operands,...
                                                     thisExp,obj.fastRedundancyCheck);
                else
                    operands=full([instrXp(row,col),instrD(col),reshape([instrLDL(row,is);instrV(is)],1,[])]);
                    instrLDL(row,col)=newInstruction(obj,obj.Itypes.I_plus_minus_dot_div,[],operands,...
                                                     thisExp,obj.fastRedundancyCheck);
                end
            else
                operands=full([instrD(col),reshape([instrLDL(row,is);instrV(is)],1,[])]);
                instrLDL(row,col)=newInstruction(obj,obj.Itypes.I_minus_dot_div,[],operands,...
                                                 thisExp,obj.fastRedundancyCheck);
            end
            instrLDL(col,row)=instrLDL(row,col);
        end
    end    

    % put diagonal back in place
    instrLDL=spdiags(instrD',0,instrLDL);
    %full(instrLDL)
    
    % keep subscripts in subsLDL in the "natural" order (sorted by row and then col)
    [i,j,instrLDL]=find(instrLDL);
    [subsLDL,k]=sortrows([uint64(i),uint64(j)],2:-1:1);
    subsLDL=subsLDL';
    instrLDL=instrLDL(k);

    if verboseLevel>1
        LDL=sparse(double(subsLDL(1,:)),...
                   double(subsLDL(2,:)),...
                   instrLDL,n,n);
        k=subsLDL(1,:)>=subsLDL(2,:);
        LD=sparse(double(subsLDL(1,k)),...
                  double(subsLDL(2,k)),...
                  instrLDL(k),n,n);
        
        figure(101);
        subplot(2,3,5)
        plot(1:n,1:n,'rx');
        legend('pivots');
        hold on;
        spy(LD);
        title('LD (csparse''s ldl (pivots in red))');
        if verboseLevel>2
            fprintf('LDL (csparse''s ldl (pivots in red))\n');
            disp(full(LDL(1:10,1:10)))
        end

        subplot(2,3,6)
        plot(1:length(p),p,'rx');
        axis ij;axis image
        legend('row perm');

        % plot pivots in original matrix
        subplot(2,3,1);
        hold on;
        plot(p,p,'rx');
        legend('pivots')

        drawnow

        %fprintf('\n<PAUSED at sparsity_lu.m>\n')
        %pause;


        %fprintf('\nq=%s;\n',mat2str(q));
        %fprintf('\np=%s;\n',mat2str(p));
        try
            figure(fig);
        catch me
            fprintf('Unable to go back to current figure (did you close the window?)\n');
        end
    end
    
    % profile off
    % profile viewer
    
    obj.statistics.ldl{end+1,1}=struct(...
        'sizeA',osize,...
        'nnzA',length(instrX),...
        'uniqueNnzA',length(unique(instrX)),...
        'nNonSymA',nNonSymA,...
        'nnzLDL',length(instrLDL),...
        'uniqueNnzLDL',length(unique(instrLDL)),...
        'nIntrLDL',instructionsTableHeight()-n0);

    if verboseLevel>0
        fprintf('  sparsifyLDL(%3d): LDL size=%-10s, nnz=%d,                                         # new instr=%4d (%d..%d) (%.2f sec)\n',...
                thisExp,['[',index2str(osize),']'],length(instrLDL),...
                instructionsTableHeight()-n0,n0+1,instructionsTableHeight(),etime(clock,t0));
    end

end

function nNonSym=symmetric(obj,subsX,instrX)
    
    subsXt=subsX([2,1],:);
    [~,k1]=sortrows(subsX');
    [~,k2]=sortrows(subsXt');

    %[subsX(:,k1)',subsXt(:,k2)']
    %[instrX(k1),instrX(k2)]
    k=find(any(subsX(:,k1)~=subsXt(:,k2),1)' | instrX(k1)~=instrX(k2));
     
    if false %~isempty(k)
        [subsX(:,k1(k))',instrX(k1(k)),subsX(:,k2(k))',instrX(k2(k)),k]
        
        i=1;
        fprintf('%d:\n',instrX(k1(k(i))));
        showInstruction(obj,instrX(k1(k(i))),2)
        fprintf('%d:\n',instrX(k2(k(i))));
        showInstruction(obj,instrX(k2(k(i))),2)
    end
    
    nNonSym=length(k);
end

function test()

    n=6;
    A=randn(n,n);
    A=A*diag(randn(n,1))*A';
    
    LDL=nan(n);
    D=nan(n,1);
    
    % without pivoting
    for col=1:n 
        v=nan(col-1,1);
        for i=1:col-1 
            v(i)=LDL(col,i)*D(i);
        end 
        D(col)=A(col,col)- LDL(col,1:col-1)*v(1:col-1);
        if D(col)==0
            error('needs pivoting');
        end

        for row=col+1:n
            LDL(row,col)=( A(row,col)- LDL(row,1:col-1)*v(1:col-1) ) / D(col) ;
            LDL(col,row)=LDL(row,col); % not needed, but just to get L'
        end
        
    end

    LDL
    L=tril(LDL,-1)+eye(n)
    D=diag(D)

    L*D*L'-A
    
    
    
end

