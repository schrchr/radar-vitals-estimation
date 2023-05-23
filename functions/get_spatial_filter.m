function [b] = get_spatial_filter(a,d_p,theta_p)
%Get_SPATIAL_FILTER outputs the spatial filters (MVDR) for each given 
% location. The filter 
%
%   Detailed explanation goes here

P = length(d_p);

if P ~= length(theta_p)
    error('d_p has not the same size as theta_p')
end

a_dummy = a(1,1);
b = NaN(length(a_dummy(:)),P);
for p = 1:P
    steer_p = a(d_p(p),theta_p(p));
    b(:,p) = steer_p(:);
end

end