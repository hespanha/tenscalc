function str=index2str(index,format,compress)
% str=index2str(index)
%   Converts an array to string
%
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


    if nargin<2
        format='%d';
    end

    if nargin<3
        compress=false;
    end

    if size(index,1)==1
        di=diff(index);
        if compress && length(index)>3 && all(di==di(1))
            str=num2str([index(1),di(1),index(end)],[format,':',format,':',format]);
        else
            str=num2str(index,[format,',']);
            if ~isempty(str)
                str(end)=[];
            end
            str=regexprep(str,' ','');
        end
    else
        str=cellstr(num2str(index,[format,',']))';
        str=strjoin(str,';');
        str=regexprep(str,',;',';');
        str=regexprep(str,',$| ','');
    end

    % if length(index)<=1
    %     str=sprintf(format,index);
    % else
    %     str=sprintf(['%s',format],sprintf([format,','],index(1:end-1)), ...
    %                 index(end));
    % end
end
