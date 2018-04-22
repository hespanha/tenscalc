function str=mat2str_compact(mat)
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
    
    msize=size(mat);

    if nnz(mat)==0
        %mat=full(mat);
        str=sprintf('zeros(%s)',index2str(msize));
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
