function [K] = RD(lambda, alpha, plt)
% Model order estimation based on the relative distance (RD) between eigenvalues
% https://ieeexplore.ieee.org/document/6145674
%
% Input:
%   lambda: (Nx1) vector of eigenvalues sorted in descending order (highest to lowest)
%   alpha: (1x1) tuning parameter (paper suggests between 2 and 5)
%   plt: (bool) turn on/off plot
%
% Output:
%   K: (1x1) estimated model order

if nargin == 1
    alpha = 3;
    plt = false;
elseif nargin == 2
    plt = false;
end

% calculate relative distance
RD_value = abs(diff(lambda))./lambda(2:end);
% get index of 5 highest RDs
[M, I] = maxk(RD_value, 5);

[I_sort,II] = sort(I,1,"descend");
M_sort = M(II);

ii = 1;
while lambda(I_sort(ii)) < alpha * mean(lambda(I_sort(ii)+1:end))
    ii = ii + 1;
    if ii > 5
        K = nan;
        warning("RDI value larger than 5, K = -1;")
        return;  
    end
end

K = round(I_sort(ii)/2);

if plt == true
    figure
    plot(RD_value)
    hold on
    plot(I_sort,M_sort,'rx')
    plot(I_sort(ii),M_sort(ii),'go')
    grid on
    xlabel('RDI')
    ylabel('RD')
    %xlim([0 50])
end

end

