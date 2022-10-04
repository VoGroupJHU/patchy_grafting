function pts_full = gen_poly(pts_core,grid_size)
% n_perm = 1E3;
% for i = 1:n_perm
%     indx_perm = randperm(length(pts_core(:,1)));
%     pts_core = pts_core(indx_perm,:);
% end
% Unique points
P = unique(pts_core,'rows');

% Generate convex hull
k = boundary(P,0);

% Points generator
a1 = linspace(0,1,grid_size+2);
a2 = linspace(0,1,grid_size+2);
% Triangulation to fill in points
pts_full = [];
for kk = 1:length(k(:,1))
    triangle_test = k(kk,:);
    polygon = P(triangle_test,:);
    % Shift to origin
    polygon_shift = zeros(length(polygon(:,1)),3);
    for i = 1:length(polygon(:,1))
        polygon_shift(i,:) = polygon(i,:) - polygon(2,:);
    end
    % Leg vectors
    v1 = polygon_shift(1,:);
    v2 = polygon_shift(3,:);
    % Generate grid
    pts_in = zeros(length(a1)*length(a2),3);
    counter = 1;
    for i = 1:length(a1)
        for j = 1:length(a2)
            p_temp = a1(i)*v1 + (1-a1(i))*a2(j)*v2;
            pts_in(counter,:) = p_temp + polygon(2,:);
            counter = counter + 1;
        end
    end
    % Combine and store
    pts_full = [pts_full; pts_in];
end