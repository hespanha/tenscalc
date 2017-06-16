function writeCasmLBinstructions(obj,fid,ks)
% Type of code produced:'
%    ''C+asmLB'' - little C code, with most of the computations done'
%            by large blocks of inlined assembly code'
%            Impact on non-optimized compilation (ideal for testing):'
%              + fastest compilation times'
%              + smallest code size'
%              + fastest run times'
%            Impact on optimized compilation:'
%              Most of the compiler optimization is restricted to re-ordering'
%              and/or inlining the large blocks of asm code'
%              . fastest compile optimization times'
%              . slowest run times'
%              . largest optimized code sizes (due to inlining large blocks)'
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

constants=[];
for i=1:length(ks)
    k=ks(i);
    [type,parameters,operands]=getInstruction(int64(k));
    if type==obj.Itypes.I_load %|| type==obj.Itypes.I_set
        constants=[constants,parameters];
    end
end

if isempty(constants)
    fprintf(fid,'\t{\n',mymat2str(constants));
else
    fprintf(fid,'\t{ double c[]={%s};\n',mymat2str(constants));
end
fprintf(fid,'\tasm(\n');
ic=0;

for i=1:length(ks)
    k=ks(i);
    [type,parameters,operands]=getInstruction(int64(k)); % 0-based operands
    operands=obj.memoryLocations(operands)-1; % -1 for 0-base indexing
    switch (type)
      case obj.Itypes.I_plus_minus_dot
        % C
        % fprintf(fid,'\t// m[%d]=m[%d]-(',obj.memoryLocations(k-1),SIZE*operands(1));
        % fprintf(fid,'+m[%d]*m[%d]',SIZE*operands(2:end));
        % fprintf(fid,'); //plus-dot\n');
        % ASM
        fprintf(fid,'\t"movsd %d(%%[m]),%%%%xmm0;',SIZE*operands(1));
        fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %d(%%[m]),%%%%xmm1;subsd %%%%xmm1,%%%%xmm0;',SIZE*operands(2:end));
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);" //plus-dot\n',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_sum
        % C
        % fprintf(fid,'\t// m[%d]=',obj.memoryLocations(k-1));
        % par=sprintf('%+dm[%d]',[parameters;SIZE*double(operands)]); % double needed since operands are uint64
        % par=regexprep(par,'([+-])1','$1');
        % fprintf(fid,'%s',par);
        % fprintf(fid,'; //sum\n');
        % ASM
        fprintf(fid,'\t"');
        %% Group summations and subtractions -- NOT GOOD SINCE AMPLIFIES NUMERICAL ERRORS
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
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);" //sum\n',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_sumprod
        % C
        % str=regexprep(index2str(SIZE*operands','m[%d]'),',','*');
        % str=regexprep(str,';','+');
        % fprintf(fid,'\t// m[%d]=%s; //sumprod\n',obj.memoryLocations(k-1),str);
        % ASM
        nSum=parameters(2);
        nProd=parameters(1);
        operands=reshape(operands,nProd,nSum);
        fprintf(fid,'\t"movsd %d(%%[m]),%%%%xmm0;',SIZE*operands(1,1));
        if nProd>1
            fprintf(fid,'mulsd %d(%%[m]),%%%%xmm0; ',SIZE*operands(2:end,1));
        end
        if nSum>1
            format=['movsd %d(%%[m]),%%%%xmm1;',...
                    repmat('mulsd %d(%%[m]),%%%%xmm1;',1,nProd-1),...
                    'addsd %%%%xmm1,%%%%xmm0; '];
            fprintf(fid,format,SIZE*operands(:,2:end));
        end
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);" //sumprod\n',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_div
        % C
        % fprintf(fid,'\t// m[%d]=m[%d]/m[%d]; //div\n',obj.memoryLocations(k-1),SIZE*operands(1),SIZE*operands(2));
        % ASM
        fprintf(fid,'\t"movsd %d(%%[m]),%%%%xmm0;divsd %d(%%[m]),%%%%xmm0;',SIZE*operands(1),SIZE*operands(2));
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);" //div\n',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_minus_dot
        % C
        % fprintf(fid,'\t// m[%d]=-(',obj.memoryLocations(k-1));
        % fprintf(fid,'m[%d]*m[%d]',operands);
        % fprintf(fid,'); //-dot\n');
        % ASM
        fprintf(fid,'\t"subsd %%%%xmm0,%%%%xmm0;');
        fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %d(%%[m]),%%%%xmm1;subsd %%%%xmm1,%%%%xmm0;',SIZE*operands);
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);" //-dot\n',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_set
        % Do nothing, since set function takes care of that
      case obj.Itypes.I_load
        % C
        % fprintf(fid,'\t// m[%d]=%g; //load\n',SIZE*obj.memoryLocations(k-1),parameters);
        % ASM
        fprintf(fid,'\t"movsd %d(%%[c]),%%%%xmm0;movsd %%%%xmm0,%d(%%[m]);" //load\n',SIZE*ic,SIZE*obj.memoryLocations(k-1));
        ic=ic+1;
      case obj.Itypes.I_minus_dot_div
        % C
        % fprintf(fid,'\tm[%d]=-(',SIZE*obj.memoryLocations(k-1));
        % fprintf(fid,'m[%d]*m[%d]',SIZE*operands(2:end));
        % fprintf(fid,')/m[%d]; //-dot-div\n',SIZE*operands(1));
        % ASM
        fprintf(fid,'\t"movsd %d(%%[m]),%%%%xmm0;mulsd %d(%%[m]),%%%%xmm0;',SIZE*operands(2:3));
        if length(operands)>3
            fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %d(%%[m]),%%%%xmm1;addsd %%%%xmm1,%%%%xmm0;',SIZE*operands(4:end));
        end
        fprintf(fid,'movsd %%%%xmm1,%%%%xmm1;subsd %%%%xmm0,%%%%xmm1;divsd %d(%%[m]),%%%%xmm1;',SIZE*operands(1));
        fprintf(fid,'movsd %%%%xmm1,%d(%%[m]);" //-dot-div\n"',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_plus_sqr
        % C
        % fprintf(fid,'\t// m[%d]=',obj.memoryLocations(k-1));
        % fprintf(fid,'+m[%d]*m[%d]',repmat(SIZE*operands(:)',2,1));
        % fprintf(fid,'; //plus-sqr\n');
        % ASM
        fprintf(fid,'\t"movsd %d(%%[m]),%%%%xmm0;mulsd %%%%xmm0,%%%%xmm0;',SIZE*operands(1));
        fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %%%%xmm1,%%%%xmm1;addsd %%%%xmm1,%%%%xmm0;',SIZE*operands(2:end));
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);" //plus-sqr\n',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_plus_minus_dot_div
        % C
        % fprintf(fid,'\t// m[%d]=(m[%d]-(',obj.memoryLocations(k-1),SIZE*operands(1));
        % fprintf(fid,'+m[%d]*m[%d]',SIZE*operands(3:end));
        % fprintf(fid,'))/m[%d]; //plus-dot-div\n',SIZE*operands(2));
        % ASM
        fprintf(fid,'\t"movsd %d(%%[m]),%%%%xmm0;mulsd %d(%%[m]),%%%%xmm0;',SIZE*operands(3:4));
        if length(operands)>4
            fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;mulsd %d(%%[m]),%%%%xmm1;addsd %%%%xmm1,%%%%xmm0;',SIZE*operands(5:end));
        end
        fprintf(fid,'movsd %d(%%[m]),%%%%xmm1;subsd %%%%xmm0,%%%%xmm1;divsd %d(%%[m]),%%%%xmm1;',SIZE*operands(1),SIZE*operands(2));
        fprintf(fid,'movsd %%%%xmm1,%d(%%[m]);" //plus-dot-div\n',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_min
        % C
        % fprintf(fid,'\t// m[%d]=m[%d];\n',SIZE*obj.memoryLocations(k-1),SIZE*operands(1));
        % if length(operands)>1
        %    fprintf(fid,'\t// if (m[%d]>m[%d]) m[%d]=m[%d]; //min\n',...
        %            [repmat(k-1,1,length(operands)-1);...
        %             SIZE*operands(2:end)';...
        %             repmat(k-1,1,length(operands)-1);...
        %             SIZE*operands(2:end)']);
	% end
        % ASM
        fprintf(fid,'\t"movsd %d(%%[m]),%%%%xmm0;',SIZE*operands(1));
        if length(operands)>1
            fprintf(fid,'minsd %d(%%[m]),%%%%xmm0;',SIZE*operands(2:end));
        end
        fprintf(fid,'movsd %%%%xmm0,%d(%%[m]);" //min\n',SIZE*obj.memoryLocations(k-1));
      case obj.Itypes.I_clp
        % C
        % fprintf(fid,'\t// m[%d]=DBL_MAX;\n',SIZE*obj.memoryLocations(k-1));
        % if size(operands,2)>0
        %     fprintf(fid,'\t// if (m[%d]<0) { typeof(*m) x=-m[%d]/m[%d]; if (m[%d]>x) m[%d]=x; } //clp\n',...
        %             [SIZE*operands(2,:);SIZE*operands(1,:);SIZE*operands(2,:);...
        %              repmat(SIZE*obj.memoryLocations(k-1),2,length(operands))]);
        % end
        % ASM
        operands=reshape(operands,2,length(operands)/2);
        fprintf(fid,'\t"movsd %%[dbl_max],%%%%xmm0;\\\n');
        fprintf(fid,'\t\t%d%d:movsd %d(%%[m]),%%%%xmm1;subsd %%%%xmm2,%%%%xmm2;ucomisd %%%%xmm2,%%%%xmm1;jae %d%df;subsd %%%%xmm1,%%%%xmm2;movsd %d(%%[m]),%%%%xmm1;divsd %%%%xmm2,%%%%xmm1;minsd %%%%xmm1,%%%%xmm0;\\\n',...
                [SIZE*obj.memoryLocations(k-1)*ones(1,size(operands,2));
                 0:size(operands,2)-1;
                 SIZE*operands(2,:);
                 SIZE*obj.memoryLocations(k-1)*ones(1,size(operands,2));
                 1:size(operands,2);
                 SIZE*operands(1,:)]);
        fprintf(fid,'\t%d%d:movsd %%%%xmm0,%d(%%[m]);" //clp\n',SIZE*obj.memoryLocations(k-1),size(operands,2),SIZE*obj.memoryLocations(k-1));
      otherwise
        obj.Itypes
        error('instruction %d not implemented\n',type)
    end
    
end
if isempty(constants)
    fprintf(fid,'\t::[m] "r" (scratchbook), [dbl_max] "m" (dbl_max) :"xmm0","xmm1","xmm2");\n');
else
    fprintf(fid,'\t::[m] "r" (scratchbook), [dbl_max] "m" (dbl_max), [c] "r" (c) :"xmm0","xmm1","xmm2");\n');
end
fprintf(fid,'\t};\n');

