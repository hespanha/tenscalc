function varargout=precompute(varargin)
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    maxNumel2ExpandConstants=2000;
    maxNnz2ExpandConstants=1000;

    varargout=varargin;

    for i=1:length(varargin)
        type=getType(varargin{i});
        if ~isequal(type,'constant') && ~isequal(type,'zeros') && ~isequal(type,'ones') && getIsStatic(varargin{i})
            osize=size(varargin{i});
            if prod(osize)>maxNumel2ExpandConstants
                fprintf('precompute: not pre-computing constant [%s] with two many elements (numel = %d>%d)\n',index2str(osize),prod(osize),maxNumel2ExpandConstants);
                continue;
            end
            x=optimize4matlab(varargin{i});
            ex=eval(str(x));
            if 0%length(size(varargin{i}))>2
                fprintf('precompute: not pre-computing constant [%s] with dimension>=2 (numel = %d, nnz = %d)\n',index2str(size(x)),prod(size(x)),nnz(ex));
                continue;
            end
            if nnz(ex)>maxNnz2ExpandConstants
                fprintf('precompute: not precomputing constant [%s] with too many non-zero elements (numel = %d, nnz = %d>%d)\n',index2str(size(x)),prod(size(x)),nnz(ex),maxNnz2ExpandConstants);
                continue;
            end
            %varargin{i}
            varargout{i}=Tconstant(ex,osize);
            fprintf('precompute: precomputing (osize: %s -> %s, nnz: %d)\n',index2str(osize),index2str(size(varargout{i})),nnz(ex));
        end
    end
end
