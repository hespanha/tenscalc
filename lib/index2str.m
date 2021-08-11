function str=index2str(index,format,compress)
% str=index2str(index)
%   Converts an array to string
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

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
