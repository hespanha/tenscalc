function br=bitrate(snr)
% br=bitrate(snr)
%
% Returns Shannon's error-free bit-rate per Hz, for an analog
% communication channel subject to additive white Gaussian noise with
% the given signal-to-noise ration: 
%     br=log2(1+snr)

br=compose(snr,@(x__)log(1+x__)/log(2),@(x__)1./(1+x__)/log(2),@(x__)-1./(1+x__).^2);
