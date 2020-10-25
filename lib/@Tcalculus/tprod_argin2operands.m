function [tprod_size,sums_size,objs,inds]=tprod_argin2operands(varargin)
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

    if mod(nargin,2)
        error('tprod: number of arguments must be even\n');
    end

    tprod_size=[];
    sums_size=[];
    objs=cell(nargin/2,1);
    inds=cell(nargin/2,1);

    for i=1:2:nargin
        objs{(i-1)/2+1}=toCalculus(varargin{i});
        ind=varargin{i+1};

        % if any(diff(sort(ind))==0)
        %     ind
        %     warning('currently tprod does not support repeated indices: ind(%d)=[%s]\n',...
        %             (i+1)/2,index2str(ind));
        % end

        if ~isnumeric(ind)
            ind
            error('tprod: argument %d must be numeric\n',i+1);
        end
        if isempty(ind)
            % make empty indices uniform to make sure (my)isequal
            % works when compaaring the indices
            inds{(i-1)/2+1}=[];
        else
            inds{(i-1)/2+1}=ind;
        end
        osize=size(objs{(i-1)/2+1});

        if length(ind)~=length(osize)
            osize
            ind
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
                    objs{(i-1)/2+1},tprod_size
                    error(['incompatible sizes found in object ' ...
                           '%d, dimension %d (%d~=%d)\n'],(i+1)/2,j, ...
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
                    objs{(i-1)/2+1},sums_size
                    error(['incompatible sizes found in object ' ...
                           '%d, dimension %d (%d~=%d)\n'],(i+1)/2,j, ...
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
