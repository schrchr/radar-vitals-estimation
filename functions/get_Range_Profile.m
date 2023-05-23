function [range_profile, range_vec] = get_Range_Profile(S,deltaf,N)
% Returns range profile 
%   S (K x ...): signal matrix with K fast time samples in the first dimension
%   deltaf: frequency step size B/K
%   N: FFT size

    c = physconst('LightSpeed');
    range_profile = abs(ifft(S,N,1)); 
    range_vec = (0:N-1).'*c./(2*deltaf*N);
end