function br=bitrate(snr)
% br=bitrate(snr)
%
% Returns Shannon's error-free bit-rate per Hz, for an analog
% communication channel subject to additive white Gaussian noise with
% the given signal-to-noise ration:
%     br=log2(1+snr)
%
% This file is part of Tencalc.
%
% Copyright (C) 2012-21 The Regents of the University of California
% (author: Dr. Joao Hespanha).  All rights reserved.

if isequal(class(snr),'Tcalculus')
    br=compose(snr,@(x__)log(1+x__)/log(2),@(x__)1./(1+x__)/log(2),@(x__)-1./(1+x__).^2);
else
    br=log(1+snr)/log(2);
end
end

