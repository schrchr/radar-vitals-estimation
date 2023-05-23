function [theta, radius] = get_Ref_Circles(ref_doa, ref_range, circle_radius, ref_P)
% data for circles around reference positions
%
%   ref_doa in grad
%   radius in meters

[xc,yc] = pol2cart(ref_doa*pi/180,ref_range);
th = linspace(0,2*pi,50);
[xc2,yc2] = pol2cart(th, circle_radius);
[theta, radius] = cart2pol(repmat(xc2.',ref_P) + xc, repmat(yc2.',ref_P) + yc );

%     figure
%     polarplot(theta, radius);
%     figure
%     plot(radius, theta*180/pi,'LineWidth',2);