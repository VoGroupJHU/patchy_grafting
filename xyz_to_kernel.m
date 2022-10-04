function [theta,phi,r] = xyz_to_kernel(pts)
% Input 
% pts: points to be converted

% Output
% coordiantes in kernel space

x = pts(:,1);
y = pts(:,2);
z = pts(:,3);

theta = atan2(y,x);
phi = atan2(z,sqrt(x.^2 + y.^2));
r = sqrt(x.^2 + y.^2 + z.^2);

end 