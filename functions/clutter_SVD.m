function [Xvs, Xdc, Xnoise] = clutter_SVD(X, k1, k2)
%PREFILTER 
%   X:  Radar data KxMxL
%   k1: Number of DC eigenvectors
%   k2: Number of signal eigenvectors


Xdc = zeros(size(X));
Xvs = zeros(size(X));
Xnoise = zeros(size(X));

X = permute(X, [1,3,2]);

for m = 1:size(X,3)
    [U,S,V] = svd(X(:,:,m),"econ");
    Xdc(:,m,:)     = U(:,1:k1)*S(1:k1,1:k1)*(V(:,1:k1).');
    Xvs(:,m,:)      = U(:,k1+1:k2)*S(k1+1:k2,k1+1:k2)*(V(:,k1+1:k2).');
    Xnoise(:,m,:)   = U(:,k2+1:end)*S(k2+1:end,k2+1:end)*(V(:,k2+1:end).');
end

end
