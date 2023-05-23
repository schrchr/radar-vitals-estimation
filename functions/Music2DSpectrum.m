function [P_MUSIC] = Music2DSpectrum(V, P, MapConstraints, a)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

delta_d = MapConstraints(1);
delta_theta = MapConstraints(2);
max_d = MapConstraints(3);
max_theta = MapConstraints(4);

Vn = V(:,P+1:end);
Vs = V(:,1:P);

Vnn = Vn * Vn';

d = 0:delta_d:max_d;
theta = -max_theta:delta_theta:max_theta;
P_MUSIC = zeros(length(d),length(theta));

for r_id = 1:length(d)
    for theta_id = 1:length(theta)
        A = a(d(r_id),theta(theta_id));
        a_steer = A(:);
        P_MUSIC(r_id,theta_id) = abs(inv(a_steer'*Vnn*a_steer));
    end
end
end