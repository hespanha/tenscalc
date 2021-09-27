function str=mat2str_compact(mat)
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

    msize=size(mat);

    if nnz(mat)==0
        %mat=full(mat);
        str=sprintf('zeros(%s)',index2str(msize));
    elseif all(mat==1,'all')
        str=sprintf('ones(%s)',index2str(msize));
    elseif length(msize)==2 & msize(2)==1 & mat(2:end)==mat(1:end-1)+1
        mat=full(mat);
        str=sprintf('((%d:%d)'')',mat(1),mat(end));
    elseif length(msize)==2 & msize(1)==1 & mat(2:end)==mat(1:end-1)+1
        mat=full(mat);
        str=sprintf('(%d:%d)',mat(1),mat(end));
    elseif length(msize)==2 & msize(2)==1 & mat(2:end)==mat(1:end-1)-1
        mat=full(mat);
        str=sprintf('((%d:-1:%d)'')',mat(1),mat(end));
    elseif length(msize)==2 & msize(1)==1 & mat(2:end)==mat(1:end-1)-1
        mat=full(mat);
        str=sprintf('(%d:-1:%d)',mat(1),mat(end));
    else
        if length(size(mat))<=2
            str=mat2str(mat);
        else
            str=serialize(mat);
            % remove final ';'
            str(end)=[];
        end
    end
end
