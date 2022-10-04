function [pos_theo,vel_gradient,acc_gradient,theta_space_poly,phi_space_poly,Fx,Fy,Fz] = parameterize_shape(pts_poly,grid_step)
% Convert to spherical
% theta_poly = atan2(pts_poly(:,2),pts_poly(:,1));
[theta_poly,phi_poly] = xyz_to_kernel(pts_poly);
% Make space
d_grid = grid_step;
shift_value = 1E-3;

% theta_space_poly = min(theta_poly)+shift_value:d_grid:max(theta_poly)-shift_value;
% phi_space_poly = min(phi_poly)+shift_value:d_grid:max(phi_poly)-shift_value;
% phi_space_poly = -pi/2:d_grid:pi/2;
% theta_space_poly = -pi:d_grid:pi;

phi_space_poly = linspace(-pi/2,pi/2,grid_step);
theta_space_poly = linspace(-pi,pi,2*length(phi_space_poly));

[theta_grid,phi_grid] = meshgrid(theta_space_poly,phi_space_poly);
% Map surface
method = 'natural';
extrap = 'nearest';
Fx = scatteredInterpolant(theta_poly,phi_poly,pts_poly(:,1),method,extrap);
Fy = scatteredInterpolant(theta_poly,phi_poly,pts_poly(:,2),method,extrap);
Fz = scatteredInterpolant(theta_poly,phi_poly,pts_poly(:,3),method,extrap);
% Generate surface
x_grid = Fx(theta_grid,phi_grid);
y_grid = Fy(theta_grid,phi_grid);
z_grid = Fz(theta_grid,phi_grid);
% Store
pos_theo = {x_grid,y_grid,z_grid};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test of gradient for superellipsoid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity gradient
vel_gradient = cell(1,length(pos_theo));
for i = 1:length(pos_theo)
    % Calculate gradient in x, y, or z direciton.99
    vel_tmp = cell(1,2);
    [vel_tmp{2},vel_tmp{1}] = gradient(pos_theo{i},d_grid,d_grid);
    % Store
    vel_gradient{i} = vel_tmp;
end
% Acceleration gradient
acc_gradient = cell(1,length(pos_theo));
for i = 1:length(vel_gradient)
    % Calculate gradient
    acc_tmp = cell(1,3);
    % d^2/dtheta^2
    acc_tmp{2} = gradient(vel_gradient{i}{2},d_grid,d_grid);
    % d^2/dphi^2
    [~,acc_tmp{1}] = gradient(vel_gradient{i}{1},d_grid,d_grid);
    % d^2/dphi/dtheta
    acc_tmp{3} = gradient(vel_gradient{i}{1},d_grid,d_grid);
    % Store
    acc_gradient{i} = acc_tmp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
