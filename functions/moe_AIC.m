function [usernumber] = moe_AIC(eig_vec,N)
%
% Function to find the number of users
% using Akaike equation:
%
% Find m, between 1 to M-1, which minimizes the
% following quantity:
%
% AIC(m) = -N(M-m) * log( g(m) / a(m) ) + m * (2M - m)
%
% g(m) is the geometric mean of the (M-m) smallest
% eigenvalues of the covariance matrix of observation. - "geomean"
% a(m) is the arithmetic mean. - function "mean"
%
[temp,M] = size(eig_vec);
M = max(temp,M);
% Putting in the ascendent order
% It only works in the ascendent form
eig_vec = sort(eig_vec,'ascend');
for ii = 1:M
    kk = M - ii + 1;
    trecho = eig_vec(1:kk);
    gii = geomean(trecho);
    aii = mean(trecho);
    % mm must vary from 0 to M - 1
    % ii starts with 1 and ends in M
    % wich means M-1 until 0
    mm = M - kk;
    AIC_values(ii) = -N * ( M - mm ) * log(gii/aii) + mm * (2*M - mm);
end
[AIC_min,position] = min( AIC_values(1:M) );
usernumber = position - 1;