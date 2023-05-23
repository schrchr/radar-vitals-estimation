function [R, W_K, W_M, W_KM] = SpacialSmoothing(S, W)
% applies spacial smoothing on S in order to reduce wave correlation
%
%   Input:
%       S   (K x M) Signal matrix
%       W   If W is scalar, the smoothing window will have approximately the same ratio W_K/W_M as S. 
%           If W is a vector, W_K and W_M can be choosen freely
%           W = [W_K, W_M].
%   Output:
%       Rxx Decorrelated covariance matrix of size (W_K W_M x W_K W_M)

    [K, M] = size(S);
 
    if max(size(W))>1
        W_K = W(1);
        W_M = W(2);
        if max(size(W))>2
            warning("Too many entries for W")
        end
        if W_K<W_M
            warning("W_K<W_M")
        end
        if W_K > K || W_M > M
            error("W_K > K || W_M > M")
        end
    else
        W_M = floor(sqrt(M/K*W));
        W_K = floor(sqrt(K/M*W));
        if W > K*M
            error("W > K*M")
        end
    end
    
    W_KM = (K - W_K + 1)*(M - W_M + 1);
    
    X = zeros(W_M * W_K, W_KM);
    
    ii = 1;
    for i = 1:(K - W_K + 1)
        for j = 1:(M - W_M + 1)
            S_ij = S(i:i+W_K-1, j:j+W_M-1); % W_K x W_M
            X(:,ii) = S_ij(:); % W_M*W_K x W_KM
            ii = ii + 1;
        end
    end
    R = X * X' ./ W_KM;
end