function [subsLU,instrLU,p,q]=sparsity_lu(obj,thisExp,typical_subscripts,typical_values)
%   Computes the sparsity pattern for an elementary expression
%   'thisExp' of type 'lu' and returns the positions in memory
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
        error('LU decomposition only implemented for square matrices\n');
    end
    
    [uniqueInstrX,~,kuniqueInstrX]=unique(instrX); % instrX=uniqueInstrX(kuniqueInstrX);

    if verboseLevel>0
        t0=clock();
        n0=instructionsTableHeight();
        fprintf('  sparsifyLU(%3d): A  size=%-10s, nnz=%d, distinct nnz=%d\n',...
                thisExp,['[',index2str(osize),']'],length(instrX),length(uniqueInstrX));
    end

    
    % compute matrix to optimize sparsity in lu factorization
    [sz,subscripts,values,A]=loadCSparse(typical_subscripts,typical_values);

    if any(isnan(sz)) || ~isequal(subscripts,subsX) 
        fprintf('    using random values, ');

        s = RandStream('mt19937ar','Seed',0);
        % create random values preserving structural symmetries
        values=1+.5*rand(s,length(uniqueInstrX),1);
        values=values(kuniqueInstrX);
        
        A=sparse(double(subsX(1,:)),double(subsX(2,:)),values,n,n);
    else
        fprintf('    using values from "%s"\n    ',typical_values);
        nnan=sum(isnan(A(:)));
        if nnan>0
            error('typical values include %d nan entries',full(nnan));
        end
        
        % add some noise to remove non structural zeros (in A its the LU factorization)
        symm=isequal(A,A');
        tol=min(find(abs(A)))*eps^2;
        A=A+sparse(double(subscripts(1,:)),double(subscripts(2,:)),tol*randn(size(subscripts,2),1));
        if symm
            % preserve symmetry
            A=(A+A')/2;
        end
    end


    nNonSymA=symmetric(obj,subsX,instrX);
    if nNonSymA==0
        fprintf('structurally SYMMETRIC matrix, ');
        if ~isequal(A,A')
            norm(full(A-A'))
            error('Unexpected non symmetric matrix');
        end
    else
        fprintf('%d/%d nonzero entries of A and A'' do not match',nNonSymA,length(instrX));
    end
    
    if 1
        % determine column permutation and pivoting for "optimal sparsity"
        % p = row permutation; q = column permutation
        [lp,up,p,q]=lu(A,.01,'vector');  % .01 gives more flexibility in minimizing fillin
    else
        
        % compute sparsifying column permutation using colamd
        q=colamd(A);
        p=q;
        %p=1:length(q);
        [lp,up]=lu(A(p,q));
    end
    
    if verboseLevel>1 && nnz(A)<prod(size(A))
        fig=get(0,'CurrentFigure');
        f=figure(101);
        set(f,'Name',typical_subscripts);
        clf

        subplot(2,4,1)
        spy(A);
        % plot(subsX(1,kSet),subsX(2,kSet),'r+')
        % plot(subsX(1,kLoad),subsX(2,kLoad),'g+')
        if nNonSymA==0
            title(sprintf('symmetric A %dx%d (original)',size(A)));
        else
            title(sprintf('asymetric A %dx%d (original)',size(A)));
        end
        hold on;
        
        if n<100
            subplot(2,4,2);
            [l,u]=lu(A);
            spy(abs(l)+abs(u));
            title('L+U (matlab''s lu without perm)')
        end

        subplot(2,4,3)
        spy(abs(lp)+abs(up))
        title('L+U (matlab''s lu with perm)')

        subplot(2,4,4)
        if isequal(p,q)
            plot(1:length(p),p,'r.');
            legend('row/col  perm','location','southoutside');
        else
            plot(1:length(p),p,'r.',1:length(q),q,'g.');
            legend('row perm','col perm','location','southoutside');
        end
        axis ij;axis image;
        
        subplot(2,4,5)
        spy(A(p,q))
        title('A (permuted)')
        % subplot(2,4,6)
        % [l,u]=lu(A(p,q));
        % spy(abs(l)+abs(u))
        % title('L+U (matlab''s lu without further perm)')

        drawnow
    end
 
    %% Compute instructions
    
    if any(instrX==0)
        error('sparsityLU: algorithm expects no instrX=0\n');
    end
    instrLU=sparse(double(subsX(1,:)),...
                   double(subsX(2,:)),...
                   double(instrX),n,n);
    instrLU=instrLU(p,q);
    Atyp=A(p,q);
    %save 'Atyp.mat' A p q Atyp
    
    rows2process=1:n;
    pp=nan(1,length(p));
    fprintf('\n');
    for col=1:n
        % find non-zero rows for column 'col'
        [rows,~]=find(instrLU(rows2process,col));
        if isempty(rows)
            error('sparsiy_lu: no pivot in column %d\n',col);
            continues;
        end
        rows=rows2process(rows);
        % select first non-zero pivot
        [pivot,k]=min(rows);
        rows(k)=[];
        rows2process(rows2process==pivot)=[];
        pp(col)=pivot;
        % use pivot to eliminate remaining non-zero entries in that column
        %fprintf('  0: col=%d, pivot=%d, rows=[%s]\n',col,pivot,index2str(rows));
        if mod(col,100)==1
            fprintf('  %4d/%4d: pivot(%4d,%4d)=%12g (%3d rows) %4.1fsec\n',...
                    col,n,pivot,col,full(Atyp(pivot,col)),length(rows),etime(clock(),t0));
        end
        if Atyp(pivot,col)==0
            error('unexpected zero pivot\n')
        end
        for i=1:length(rows)
            row=rows(i);
            % Get 1 in U's diagonal:
            %     LU(row,col)=LU(row,col)/LU(pivot,col);
            Atyp(row,col)=Atyp(row,col)/Atyp(pivot,col);
            instrLU(row,col)=newInstruction(obj,obj.Itypes.I_div,[],...
                                            full([instrLU(row,col),instrLU(pivot,col)]),...
                                            thisExp,obj.fastRedundancyCheck);
            %fprintf('  a: col=%d,row=%d, pivot=%d, instr=%s\n',col,row,pivot,index2str(instrLU(row,col)));
            % Subtract scaled pivot row 
            [~,cols]=find(instrLU(row,col+1:n) & instrLU(pivot,col+1:n));
            colsnz=cols+col;
            [~,cols]=find(instrLU(row,col+1:n)==0 & instrLU(pivot,col+1:n));
            colsz=cols+col;
            if ~isempty(colsnz)
                %   When LU(row,col+1:n)~=0
                %     LU(row,col+1:n)=LU(row,col+1:n)-LU(row,col)*LU(pivot,col+1:n);
                Atyp(row,colsnz)=Atyp(row,colsnz)-Atyp(row,col)*Atyp(pivot,colsnz);
                %fprintf('plus-mult\n')
                % operands=full([instrLU(row,colsnz)',...
                %                repmat(instrLU(row,col),length(colsnz),1),...
                %                instrLU(pivot,colsnz)']);
                % instrLU(row,colsnz)=newInstructions(obj,obj.Itypes.I_plus_minus_dot,{[]},...
                %                                     num2cell(operands,2),...
                %                                     thisExp,obj.fastRedundancyCheck);
                for j=1:length(colsnz)
                    instrLU(row,colsnz(j))=newInstruction(obj,...
                                                          obj.Itypes.I_plus_minus_dot,[],...
                                                          full([instrLU(row,colsnz(j)),...
                                                           instrLU(row,col),...
                                                           instrLU(pivot,colsnz(j))]),...
                                                          thisExp,obj.fastRedundancyCheck);
                end
                %fprintf('  b: col=%d,row=%d, instr=%s\n',col,row,index2str(instrLU(row,colsnz)));
            end
            if ~isempty(colsz)
                %   When LU(row,col+1:n)==0
                %     LU(row,col+1:n)=-LU(row,col)*LU(col,col+1:n);
                Atyp(row,colsz)=-Atyp(row,col)*Atyp(pivot,colsz);
                %fprintf('-mult\n')
                % operands=full([repmat(instrLU(row,col),length(colsz),1),...
                %                instrLU(pivot,colsz)']);
                % instrLU(row,colsz)=newInstructions(obj,obj.Itypes.I_minus_dot,{[]},...
                %                                    num2cell(operands,2),...
                %                                    thisExp,obj.fastRedundancyCheck);
                for j=1:length(colsz)
                    instrLU(row,colsz(j))=newInstruction(obj,...
                                                          obj.Itypes.I_minus_dot,[],...
                                                          full([instrLU(row,col),...
                                                           instrLU(pivot,colsz(j))]),...
                                                          thisExp,obj.fastRedundancyCheck);
                end
                %fprintf('  c: col=%d,row=%d, instr=%s\n',col,row,index2str(instrLU(row,colsz)));
            end
        end
    end
    if verboseLevel>1
        if isequal(pp,1:length(pp))
            fprintf('no changes in pivoting from matlab''s lu, \n');
        else
            fprintf('changing %d pivots, \n',sum(pp~=1:length(pp)));
        end
    end
    %any(pp(2:end)-pp(1:end-1)~=1)
    p=p(pp);
    instrLU=instrLU(pp,:);
    %    full(instrLU)
    
    % keep subscripts in subsLU in the "natural" order (sorted by row and then col)
    [i,j,instrLU]=find(instrLU);
    [subsLU,k]=sortrows([uint64(i),uint64(j)],2:-1:1);
    subsLU=subsLU';
    instrLU=instrLU(k);

    if verboseLevel>1 && nnz(A)<prod(size(A))
        LU=sparse(double(subsLU(1,:)),...
                  double(subsLU(2,:)),...
                  ones(1,size(subsLU,2)),n,n);
        
        figure(f);
        subplot(2,4,7)
        spy(LU);
        hold on;
        plot(1:n,1:n,'r.');
        legend('nz entries','pivots','location','southoutside');
        title('L+U (csparse''s lu (pivots in red))')
        subplot(2,4,8)
        plot(1:length(p),p,'r.',1:length(q),q,'g.');
        axis ij;axis image
        legend('row perm','col perm','location','southoutside');

        % plot pivots in original matrix
        subplot(2,4,1);
        plot(p,q,'r.');
        legend('nz entries','pivots','location','southoutside')

        drawnow

        %fprintf('\n<PAUSED at sparsity_lu.m>\n')
        %pause;


        %fprintf('\nq=%s;\n',mat2str(q));
        %fprintf('\np=%s;\n',mat2str(p));
        if ~isempty(fig)
            try
                figure(fig);
            catch me
                fprintf('Unable to go back to current figure (did you close the window?)\n');
            end
        end
    end
    
    % profile off
    % profile viewer
    
    obj.statistics.lu{end+1,1}=struct(...
        'sizeA',osize,...
        'nnzA',length(instrX),...
        'uniqueNnzA',length(unique(instrX)),...
        'nNonSymA',nNonSymA,...
        'nnzLU',length(instrLU),...
        'uniqueNnzLU',length(unique(instrLU)),...
        'nIntrLU',instructionsTableHeight()-n0);

    if verboseLevel>0
        % fprintf('  sparsifyLU(%3d): A  size=%-10s, nnz=%d\n',...
        %         thisExp,['[',index2str(osize),']'],length(instrX));
        fprintf('  sparsifyLU(%3d): LU size=%-10s, nnz=%d,                                         # new instr=%4d (%d..%d) (%.2f sec)\n',...
                thisExp,['[',index2str(osize),']'],length(instrLU),...
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

    if 0
        n=6;
        X=rand(n,n);
        X=X*X';
        LU=X;
    elseif 1
        load Atyp;
        LU=Atyp;
        n=size(LU,1);
    else
        load Atyp;
    
        [lp,up,pA,qA]=lu(A,.01,'vector');
        LU=A(pA,qA);
        
        norm(Atyp-LU,'fro');
        k=find(abs(Atyp-LU)>1);
        format shorte
        disp(full([Atyp(k),LU(k)]))
    end
    
    % without pivoting
    for col=1:n
        fprintf('  %4d pivot(%4d,%4d)=%g\n',col,col,col,full(LU(col,col)));
        for row=col+1:n
            LU(row,col)=LU(row,col)/LU(col,col);                     % finished L
            LU(row,col+1:n)=LU(row,col+1:n)-LU(row,col)*LU(col,col+1:n); % will become U
        end
    end

    % with pivoting
    rows2process=1:n;
    p=1:n;
    pp=p;
    for col=1:n
        pivot=ceil(rand*length(rows2process))
        rows2process(rows2process==pivot)=[]
        pp(col)=p(pivot)
        for row=rows2process
            LU(row,col)=LU(row,col)/LU(pivot,col);                     % finished L
            LU(row,col+1:n)=LU(row,col+1:n)-LU(row,col)*LU(pivot,col+1:n); % will become U
        end
    end
    p=pp;
    LU=LU(p,:)
    
    
    [i,j,v]=find(LU);
    k=i>j;
    L=full(sparse(i(k),j(k),v(k),n,n)+speye(n));
    k=i<=j;
    U=full(sparse(i(k),j(k),v(k),n,n));

    LU
    L
    U

    X
    L*U
    
end

