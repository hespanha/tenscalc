function bool=myisequal(a,b)
% version of the matlab function isqual that deals correctly with
% empty arrays (e.g., [] == ones(1,0))
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

if isempty(a) && isempty(b)
    bool=true;
    return
end
bool=isequal(a,b);
end
