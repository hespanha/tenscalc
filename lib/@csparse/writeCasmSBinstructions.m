function writeCasmSBinstructions(obj,fid,ks)
% Type of code produced:'
%    ''C+asmSB'' - little C code, with most of the computations done'
%            by small blocks of inlined assembly code'
%            Impact on non-optimized compilation:'
%              . medium compilation times'
%              . medium code size'
%              . medium run times'
%            Impact on optimized code:'
%              Most of the compiler optimization is restricted to re-ordering'
%              and/or inlining the small blocks of asm code'
%              . medium compile optimization times,'
%              . medium run times'
%              . medium code sizes'
%
%
% Assumptions
% - gnu extended ams
% - Intel-x86
% - scratchbook is double
% - temporary register used: xmm0 (double) 
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

%% To improve
% is there a good way to clear a register to zero? (now using subsd %%%%xmm0,%%%%xmm0;)
SIZE=8;

typeAll=getMulti(obj.instructions,'type',ks);
parametersAll=getMulti(obj.instructions,'parameters',ks);
operandsAll=getMulti(obj.instructions,'operands',ks);

constants=[];
for i=1:length(ks)
    k=ks(i);
    type=typeAll{i};
    if strcmp(type,'I_load') %|| strcmp(type,'I_set')
        parameters=parametersAll{i};
        constants=[constants,parameters];
    end
end

fprintf(fid,'\t{ double c[]={%s};\n',mymat2str(constants));
ic=0;

for i=1:length(ks)
    k=ks(i);
    type=typeAll{i};
    parameters=parametersAll{i};
    operands=operandsAll{i}-1; % -1 for 0-base indexing
    switch (type)
      case 'I_set'
      case 'I_load'
        % C
        % fprintf(fid,'\t// m[%d]=%g; //load\n',SIZE*(k-1),parameters);
        % ASM
        fprintf(fid,'\tasm ("movsd %d(%%[c]),%%%%xmm0;movsd %%%%xmm0,%d(%%[m]);"',SIZE*ic,SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook),[c] "r" (c)  : "xmm0"); //load\n');
        ic=ic+1;
      case 'I_sum'
        % C
        % fprintf(fid,'\t// m[%d]=',k-1);
        % par=sprintf('%+dm[%d]',[parameters;SIZE*double(operands)]); % double needed since operands are uint64
        % par=regexprep(par,'([+-])1','$1');
        % fprintf(fid,'%s',par);
        % fprintf(fid,'; //sum\n');
        % ASM
        fprintf(fid,'\tasm ("');
        adds=find(parameters>0);
        if isempty(adds)
            fprintf(fid,'subsd %%%%xmm0,%%%%xmm0;');
        else
            fprintf(fid,'movsd %d(%%[m]),%%%%xmm0;',SIZE*operands(adds(1)));
            if length(adds)>1
                fprintf(fid,'addsd %d(%%[m]),%%%%xmm0;',SIZE*operands(adds(2:end)));
            end
        end
        subs=find(parameters<0);
        if ~isempty(subs)
            fprintf(fid,'subsd %d(%%[m]),%%%%xmm0;',SIZE*operands(subs));
        end
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0"); //sum\n');
      case 'I_sumprod'
        % C
        % str=regexprep(index2str(SIZE*operands','m[%d]'),',','*');
        % str=regexprep(str,';','+');
        % fprintf(fid,'\t// m[%d]=%s; //sumprod\n',k-1,str);
        % ASM
        nSum=parameters(2);
        nProd=parameters(1);
        operands=reshape(operands,nProd,nSum);
        fprintf(fid,'\tasm ("movsd %d(%%[m]),%%%%xmm0;',SIZE*operands(1,1));
        if nProd>1
            fprintf(fid,'mulsd %d(%%[m]),%%%%xmm0; ',SIZE*operands(2:end,1));
        end
        if nSum>1
            format=['movsd %d(%%[m]),%%%%xmm1;',...
                    repmat('mulsd %d(%%[m]),%%%%xmm1;',1,nProd-1),...
                    'addsd %%%%xmm1,%%%%xmm0; '];
            fprintf(fid,format,SIZE*operands(:,2:end));
        end
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0","xmm1"); //sumprod\n');
      case 'I_div'
        % C
        % fprintf(fid,'\t// m[%d]=m[%d]/m[%d]; //div\n',k-1,SIZE*operands(1),SIZE*operands(2));
        % ASM
        fprintf(fid,'\tasm ("movsd %d(%%[m]),%%%%xmm0;divsd %d(%%[m]),%%%%xmm0;',SIZE*operands(1),SIZE*operands(2));
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0"); //div\n');
      case 'I_minus_dot'
        % C
        % fprintf(fid,'\t// m[%d]=-(',k-1);
        % fprintf(fid,'m[%d]*m[%d]',operands);
        % fprintf(fid,'); //-dot\n');
        % ASM
        fprintf(fid,'\tasm ("subsd %%%%xmm0,%%%%xmm0;');
        fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %d(%%[m]),%%%%xmm1;subsd %%%%xmm1,%%%%xmm0;',SIZE*operands);
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0","xmm1"); //-dot\n');
      case 'I_minus_dot-div'
        % C
        % fprintf(fid,'\tm[%d]=-(',SIZE*(k-1));
        % fprintf(fid,'m[%d]*m[%d]',SIZE*operands(2:end));
        % fprintf(fid,')/m[%d]; //-dot-div\n',SIZE*operands(1));
        % ASM
        fprintf(fid,'\tasm ("movsd %d(%%[m]),%%%%xmm0;mulsd %d(%%[m]),%%%%xmm0;',SIZE*operands(2:3));
        if length(operands)>3
            fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %d(%%[m]),%%%%xmm1;addsd %%%%xmm1,%%%%xmm0;',SIZE*operands(4:end));
        end
        fprintf(fid,'movsd %%%%xmm1,%%%%xmm1;subsd %%%%xmm0,%%%%xmm1;divsd %d(%%[m]),%%%%xmm1;',SIZE*operands(1));
        fprintf(fid,'movsd %%%%xmm1,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0","xmm1"); //-dot-div\n');
      case 'I_plus_sqr'
        % C
        % fprintf(fid,'\t// m[%d]=',k-1);
        % fprintf(fid,'+m[%d]*m[%d]',repmat(SIZE*operands(:)',2,1));
        % fprintf(fid,'; //plus-sqr\n');
        % ASM
        fprintf(fid,'\tasm ("movsd %d(%%[m]),%%%%xmm0;mulsd %%%%xmm0,%%%%xmm0;',SIZE*operands(1));
        fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %%%%xmm1,%%%%xmm1;addsd %%%%xmm1,%%%%xmm0;',SIZE*operands(2:end));
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0","xmm1"); //plus-sqr\n');
      case 'I_plus_minus_dot'
        % C
        % fprintf(fid,'\t// m[%d]=m[%d]-(',k-1,SIZE*operands(1));
        % fprintf(fid,'+m[%d]*m[%d]',SIZE*operands(2:end));
        % fprintf(fid,'); //plus-dot\n');
        % ASM
        fprintf(fid,'\tasm ("movsd %d(%%[m]),%%%%xmm0;',SIZE*operands(1));
        fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %d(%%[m]),%%%%xmm1;subsd %%%%xmm1,%%%%xmm0;',SIZE*operands(2:end));
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0","xmm1"); //plus-dot\n');
      case 'I_plus_minus_dot_div'
        % C
        % fprintf(fid,'\t// m[%d]=(m[%d]-(',(k-1),SIZE*operands(1));
        % fprintf(fid,'+m[%d]*m[%d]',SIZE*operands(3:end));
        % fprintf(fid,'))/m[%d]; //plus-dot-div\n',SIZE*operands(2));
        % ASM
        fprintf(fid,'\tasm ("movsd %d(%%[m]),%%%%xmm0;mulsd %d(%%[m]),%%%%xmm0;',SIZE*operands(3:4));
        if length(operands)>4
            fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %d(%%[m]),%%%%xmm1;addsd %%%%xmm1,%%%%xmm0;',SIZE*operands(5:end));
        end
        fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;subsd %%%%xmm0,%%%%xmm1;divsd %d(%%[m]),%%%%xmm1;',SIZE*operands(1),SIZE*operands(2));
        fprintf(fid,'movsd %%%%xmm1,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0","xmm1"); //plus-dot-div\n');
      case 'I_min'
        % C
        % fprintf(fid,'\t// m[%d]=m[%d];\n',SIZE*(k-1),SIZE*operands(1));
        % if length(operands)>1
        %    fprintf(fid,'\t// if (m[%d]>m[%d]) m[%d]=m[%d]; //min\n',...
        %            [repmat(k-1,1,length(operands)-1);...
        %             SIZE*operands(2:end)';...
        %             repmat(k-1,1,length(operands)-1);...
        %             SIZE*operands(2:end)']);
	% end
        % ASM
        fprintf(fid,'\tasm ("movsd %d(%%[m]),%%%%xmm0;',SIZE*operands(1));
        if length(operands)>1
            fprintf(fid,'minsd %d(%%[m]),%%%%xmm0;',SIZE*operands(2:end));
        end
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);"',SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook) :"xmm0"); //min\n');
      case 'I_clp'
        % C
        % fprintf(fid,'\t// m[%d]=DBL_MAX;\n',SIZE*(k-1));
        % if size(operands,2)>0
        %     fprintf(fid,'\t// if (m[%d]<0) { typeof(*m) x=-m[%d]/m[%d]; if (m[%d]>x) m[%d]=x; } //clp\n',...
        %             [SIZE*operands(2,:);SIZE*operands(1,:);SIZE*operands(2,:);...
        %              repmat(SIZE*(k-1),2,length(operands))]);
        % end
        % ASM
        operands=reshape(operands,2,length(operands)/2);
        fprintf(fid,'\tasm ("movsd %%[dbl_max],%%%%xmm0;\\\n');
        fprintf(fid,'\t\t%d%d:movsd %d(%%[m]),%%%%xmm1;subsd %%%%xmm2,%%%%xmm2;ucomisd %%%%xmm2,%%%%xmm1;jae %d%df;subsd %%%%xmm1,%%%%xmm2;movsd %d(%%[m]),%%%%xmm1;divsd %%%%xmm2,%%%%xmm1;minsd %%%%xmm1,%%%%xmm0;\\\n',...
                [SIZE*(k-1)*ones(1,size(operands,2));
                 0:size(operands,2)-1;
                 SIZE*operands(2,:);
                 SIZE*(k-1)*ones(1,size(operands,2));
                 1:size(operands,2);
                 SIZE*operands(1,:)]);
        fprintf(fid,'\t%d%d:movsd %%%%xmm0,%d(%%[m]);"',SIZE*(k-1),size(operands,2),SIZE*(k-1));
        fprintf(fid,' ::[m] "r" (scratchbook), [dbl_max] "m" (dbl_max) :"xmm0","xmm1","xmm2"); //clp\n');
        % % ASM-like C
        % fprintf(fid,'\t{ double xmm0,xmm1,xmm2;xmm0=dbl_max[0];\\\n');
        % fprintf(fid,'\t\tL%d_%d:xmm1=m[%d];xmm2=xmm2-xmm2;if (xmm1>=xmm2) goto L%d_%d;xmm2=xmm2-xmm1;xmm1=m[%d];xmm1=xmm1/xmm2;xmm0=fmin(xmm0,xmm1);\\\n',...
        %         [SIZE*(k-1)*ones(1,size(operands,2));
        %          0:size(operands,2)-1;
        %          SIZE/SIZE*operands(2,:);
        %          SIZE*(k-1)*ones(1,size(operands,2));
        %          1:size(operands,2);
        %          SIZE/SIZE*operands(1,:)]);
        % fprintf(fid,'\tL%d_%d:m[%d]=xmm0;} //clp\n',SIZE*(k-1),size(operands,2),SIZE/SIZE*(k-1));
      case 'I_abs'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=fabs(m[%d]); } //abs\n',k-1,operands(1));
        % ASM
      case 'I_sqrt'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=sqrt(m[%d]); } //sqrt\n',k-1,operands(1));
        % ASM
      case 'I_Dsqrt'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=.5*pow(m[%d],-.5); } //Dsqrt\n',k-1,operands(1));
        % ASM
      case 'I_DDsqrt'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=-.25*pow(m[%d],-1.5); } //DDsqrt\n',k-1,operands(1));
        % ASM
      case 'I_inv'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=1/m[%d]; } //inv\n',k-1,operands(1));
        % ASM
      case 'I_minus_inv_sqr'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=-1/(m[%d]*m[%d]); } //-inv-sqr\n',k-1,operands(1),operands(1));
        % ASM
      case 'I_cos'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=cos(m[%d]); } //cos\n',k-1,operands(1));
        % ASM
      case 'I_minus_cos'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=-cos(m[%d]); } //-cos\n',k-1,operands(1));
        % ASM
      case 'I_sin'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=sin(m[%d]); } //sin\n',k-1,operands(1));
        % ASM
      case 'I_minus_sin'
        % C
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=-sin(m[%d]); } //-sin\n',k-1,operands(1));
        % ASM
      case 'I_log'
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=log(m[%d]); } //log\n',k-1,operands(1));
        % ASM
      case 'I_exp'
        fprintf(fid,'\t { double *m=scratchbook;\n');
        fprintf(fid,'\tm[%d]=exp(m[%d]); } //exp\n',k-1,operands(1));
      otherwise
        error('instruction ''%s'' not implemented\n',type)
    end
    
end

fprintf(fid,'\t};\n');

