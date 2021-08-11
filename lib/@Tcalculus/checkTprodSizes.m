function [tprod_size,sums_size]=checkTprodSizes(varargin)
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.
% Copyright 2012-2017 Joao Hespanha

    if mod(nargin,2)
        error('tprod: number of arguments must be even\n');
    end

    tprod_size=[];
    sums_size=[];
    for i=1:2:nargin
        obj=varargin{i};
        ind=varargin{i+1};
        if ~isnumeric(ind)
            ind
            error('tprod: argument %d must be numeric\n',i+1);
        end
        osize=size(obj);

        if length(ind)~=length(osize)
            fprintf('tprod(');
            for i=1:2:nargin
                obj_=toCalculus(varargin{i});
                ind_=varargin{i+1};
                fprintf('M[%s],[%s],',index2str(size(obj_)),index2str(ind_));
            end
            fprintf(')\n');
            obj,ind
            error(['tprod requires each indice to have the same number of ' ...
                   'entries as the size of the corresponding object\n'])
        end

        for j=1:length(ind)
            if ind(j)>0
                if length(tprod_size)<ind(j)
                    % pad with nan
                    tprod_size(end+1:ind(j))=nan;
                end
                if length(tprod_size)>=ind(j) && tprod_size(ind(j))>0 ...
                        && tprod_size(ind(j))~=osize(j)
                    obj,tprod_size
                    error(['incompatible sizes found in object ' ...
                           '%d, dimension %d (%d~=%d)\n'],i,j, ...
                          tprod_size(ind(j)),osize(j));
                end
                tprod_size(ind(j))=osize(j);
            else
                if length(sums_size)<-ind(j)
                    % pad with nan
                    sums_size(end+1:-ind(j))=nan;
                end
                if length(sums_size)>=-ind(j) && sums_size(-ind(j))>0 ...
                        && sums_size(-ind(j))~=osize(j)
                    obj,sums_size
                    error(['incompatible sizes found in object ' ...
                           '%d, dimension %d (%d~=%d)\n'],i,j, ...
                          sums_size(-ind(j)),osize(j));
                end
                sums_size(-ind(j))=osize(j);
            end
        end
    end

    if any(isnan(tprod_size)) || any(isnan(sums_size))
        tprod_size,sums_size
        error('tprod has no size for some indices/summations\n',i)
    end
end
