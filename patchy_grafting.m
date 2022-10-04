clear; close all; clc;

warning off;

%%%%%%%%%%%%%%%%%%%%%%%
%%% Define colormap %%%
%%%%%%%%%%%%%%%%%%%%%%%

n_color = 1000;
color_mtx = [0 3 169; 244 244 244; 174 0 0]/255;
n_color_type = length(color_mtx(:,1));
n_color_grid = ceil(n_color/n_color_type);
mymap = [];
for i = 1:n_color_type-1
    tmp_ith = color_mtx(i,:);
    tmp_jth = color_mtx(i+1,:);
    tmp_mtx = [linspace(tmp_ith(1),tmp_jth(1),n_color_grid)', linspace(tmp_ith(2),tmp_jth(2),n_color_grid)', linspace(tmp_ith(3),tmp_jth(3),n_color_grid)'];
    mymap = [mymap; tmp_mtx];
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Grafting parameter %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define in-sphere radius
ro = 5.0;

% Kuhn length
b = 1.5;
beta_kuhn = b;

% Excluded volume
v = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Grow particle to appropriate size %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define shape
shape = 'Prism3';
pts_coord = textread(strcat(shape,'.txt'));
pts_coord = pts_coord(:,1:3);

% Shorten prism height
if strcmp(shape,'Prism3') == 1
    side = norm(pts_coord(1,:)-pts_coord(3,:));
    height = range(pts_coord(:,3));
    exp_ratio = 21.6/58.9;
    pts_coord(:,3) = pts_coord(:,3)*exp_ratio;
end

[~,r_com] = calc_inertial_tensor(pts_coord);
pts_coord = bsxfun(@minus,pts_coord,r_com);

% Generate sphere
[xs,ys,zs] = sphere;
xs = ro*xs;
ys = ro*ys;
zs = ro*zs;

% Initialize
flag_out = 1;
pts_tmp = pts_coord;
factor = 1.01;
factor_tracking = 1;

% Grow shape to match insphere radius
while flag_out == 1
    disp(factor_tracking)
    tmp_shape = alphaShape(pts_tmp);
    indx_in = inShape(tmp_shape,xs,ys,zs);
    indx = find(indx_in == 1);
    if length(indx) == length(xs(:))
        flag_out = 0;
    else
        pts_tmp = pts_tmp*factor;
        factor_tracking = factor_tracking*factor;
    end
end
% Compute area
poly_area = surfaceArea(tmp_shape);

% Update
pts_coord = factor_tracking*pts_coord;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define core shape grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate surface
grid_size = 50;
pts_poly = gen_poly(pts_coord,grid_size);
pts_poly = roundn(pts_poly,-3);
pts_poly = unique(pts_poly,'rows');

% Reduce number of points
load('core_data.mat') % For quicker remove overlap

%%% CHANGE THESE %%%
sigma = 0.1;    % Grafting density
N = 50;         % Chain length    
kT = 1.0;       % Temperature
%%%%%%%%%%%%%%%%%%%%

% Determing number of grafts
fo = sigma;
fo_array = (1:numel(pts_overlap(:,1)))*(pi*b^2)/poly_area;
[~,indx_fo] = min(abs(fo_array-fo));
n_graft = indx_fo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create augmented shape kernel %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine kernel parameterization
grid_size = 100;
pts_poly = gen_poly(pts_coord/factor_tracking,grid_size);
pts_poly = roundn(pts_poly,-2);
pts_poly = unique(pts_poly,'rows');
[pos_theo,~,~,theta_space_poly,phi_space_poly,Fx,Fy,Fz] = parameterize_shape(pts_poly,grid_size);

% Define step size
dtheta = theta_space_poly(2)-theta_space_poly(1);
dphi = phi_space_poly(2)-phi_space_poly(1);

% Compute kernel
x_grid = pos_theo{1};
y_grid = pos_theo{2};
z_grid = pos_theo{3};
pts_grid = [x_grid(:) y_grid(:) z_grid(:)];

% Calculate Omega
omega = sqrt(x_grid.^2 + y_grid.^2 + z_grid.^2);
omega = omega/min(omega(:));

% Generate meshgrid
[tt,pp] = meshgrid(theta_space_poly,phi_space_poly);

% Replicate image array
rep_img = [-1 0 1];

% Create periodic images
kernel_img = [];
% Loop for theta
for i = 1:length(rep_img)
    % Image
    theta_img = rep_img(i);
    % Replication
    tt_rep = tt + theta_img*(2*pi);
    % Loop for phi
    for j = 1:length(rep_img)
        phi_img = rep_img(j);
        pp_rep = pp + phi_img*pi;
        % Storing
        tmp = [tt_rep(:) pp_rep(:) omega(:)];
        kernel_img = [kernel_img; tmp];
    end
end

% Remove overlap points
kernel_img = roundn(kernel_img,-4);
[~,indx_img] = unique(kernel_img(:,1:2),'stable','rows');
kernel_img = kernel_img(indx_img,:);

% Define bracket range to use
%%% NOTE: Take n steps back and forward in each direction %%%
n_step = 3;
% phi range
phi_low = -pi/2 - n_step*dphi;
phi_high = pi/2 + n_step*dphi;
% theta range
theta_low = -pi - n_step*dtheta;
theta_high = pi + n_step*dtheta;
% Find inde
indx = find( (kernel_img(:,1) >= theta_low) & (kernel_img(:,1) <= theta_high) & (kernel_img(:,2) >= phi_low) & (kernel_img(:,2) <= phi_high));

% Grab relevant point
kernel_table = kernel_img(indx,:);
kernel_table = sortrows(kernel_table, [2 1]);

% Compute omega for working overlap shape
[theta_overlap,phi_overlap] = xyz_to_kernel(pts_overlap);
kernel_overlap = griddata(kernel_table(:,1),kernel_table(:,2),kernel_table(:,3),theta_overlap,phi_overlap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-grafting process %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of Monte Carlo runs
n_mc = 100;

% Define space to loop over
N_array = N;
kT_array = kT;

% Create structure
[NN,TT] = meshgrid(N_array,kT_array);
phase_space = [NN(:) TT(:)];

% Initialize
corona_store = cell(1,numel(NN));
graft_store = cell(1,numel(NN));
counter = 0;

% Compute critical chi
sigma_chi = sigma*(4*pi*b^2);
chi_star = 0.5*(N - sigma_chi.^(1/4).*b.^(1/2).*max(kernel_table(:)).^(-3/4).*v^(-1).*(b./ro).^(-1/2).*sqrt(N))/N;
chi = chi_star;

% Define input
if chi/kT > (0.5-1E-3)      % For numerical stability
    vargin = [ro sigma b N (0.5-1E-3)];
else
    vargin = [ro sigma b N chi/kT];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%% Grafting %%%
%%%%%%%%%%%%%%%%

%%% Reassign parameters %%%
% Core size
ro = vargin(1);
% Grafting density
sigma = vargin(2);
fo = sigma;
% Kuhn length
b = vargin(3);
beta_kuhn = b;
% Excluded volume
v = 1.0;
% Chain length
N = vargin(4);
% Interaction
chi = vargin(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define core
pts_core = pts_overlap;

% Effective core
pts_poly = pts_core;

% Kernel
kernel = sqrt(sum(pts_core.^2,2));
kernel_scaled = kernel/max(kernel);
kernel = kernel/min(kernel);

% Mean number of distance
rdist_scale = mean(pdist(pts_core));

% Identify graftin point
[~,indx_0] = max(abs(kernel));

% Begin MC
% Performing grafting Monte Carlo
graft_count = zeros(numel(pts_poly(:,1)),1);
R_graft = zeros(length(pts_poly(:,1)),1);
P_mc_store = cell(1,n_mc);
E_mc_store = cell(1,n_mc);
E_entropic_mc_store = cell(1,n_mc);
for ii = 1:n_mc
    fprintf('N: %i, MC Run: %i out of %i\n',N,ii,n_mc)
    % Grafting points matrix
    graft_indx = zeros(n_graft,1);
    P_store = cell(1,n_graft);
    E_store = cell(1,n_graft);
    E_entropic_store = cell(1,n_graft);
    for j = 1:n_graft
%         fprintf('N: %i, MC Run: %i, Chain: %i\n',N,ii,j)
        flag_add = 0;
        while flag_add ~= 1
            % Grafting probability
            p_jth = rand;            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Generate grafting probability %
            % Select grafting location
            if j == 1
                omega_graft = zeros(size(kernel));
                omega_graft(indx_0) = 1;%kernel_scaled(indx_0);
            else
                % Select for grafts
                indx_tmp = find(graft_indx ~= 0);
                pts_grafts = pts_poly(graft_indx(indx_tmp),:);
                Rg_grafts = Rg_chain(graft_indx(indx_tmp));
                indx_full = [];
                % Find points within grafts
                for gg = 1:length(indx_tmp)
                    dr_tmp = bsxfun(@minus,pts_poly,pts_grafts(gg,:));
                    dr_tmp = sqrt(sum(dr_tmp.^2,2))-rdist_scale;
                    indx_tmp = find(dr_tmp < 2.0*(Rg_grafts(gg)-ro));
                    indx_full = [indx_full; indx_tmp];
                end
                indx_full = unique(indx_full);
                omega_graft = zeros(size(kernel));
                omega_graft(indx_full) = 1;%kernel_scaled(indx_full);
            end            
            % pdf handle functions
            v_chain = @(omega,omega_graft,chi) (v*sqrt((ones(size(omega))-2*chi*omega_graft).^2)).^(5/5);
            v_graft = v_chain(kernel,omega_graft,chi);
            R_chain = @(omega,N,omega_graft,chi) ro*fo^(1/5)*(v_graft).^(1/5).*beta_kuhn^(2/5)*N^(3/5).*omega.^(-3/5)*(beta_kuhn/ro)^(3/5) - 0.0*ro;
            R_chain_nocore = @(omega,N,omega_graft,chi) ro*fo^(1/5)*(v_graft).^(1/5).*beta_kuhn^(2/5)*N^(3/5).*omega.^(-3/5)*(beta_kuhn/ro)^(3/5) - 0.0*ro;
            P_chain = @(omega,N,omega_graft,chi) exp( -R_chain(omega,N,omega_graft,chi).^2/(N*beta_kuhn^2) - v_graft.*fo*ro^2*N^2./((R_chain(omega,N,omega_graft,chi)+1*ro).*omega).^3);
            E_chain = @(omega,N,omega_graft,chi) (-R_chain(omega,N,omega_graft,chi).^2/(N*beta_kuhn^2) - v*fo*ro^2*N^2./((R_chain(omega,N,omega_graft,chi)+ro).*omega).^3);            
            E_entropic = @(omega,N,omega_graft,chi) (-R_chain(omega,N,omega_graft,chi).^2/(N*beta_kuhn^2));            
            E_enthalpic = @(omega,N,omega_graft,chi) (-v*fo*ro^2*N^2./((R_chain(omega,N,omega_graft,chi)+ro).*omega).^3);            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute chain Rg
            Rg_chain = R_chain_nocore(kernel,N,omega_graft,chi)/6;
            
            % Compute grafting density
            P_graft = P_chain(kernel,N,omega_graft,chi);      
            
            % Rescale
            graft_filled = graft_indx(graft_indx > 0);
            indx_unique = unique(graft_filled);           
            P_graft = (P_graft-min(P_graft))/(max(P_graft)-min(P_graft));
            
            % Chain energy
            E_graft = E_chain(kernel,N,omega_graft,chi);
            
            % Entropic component
            %%% Zero out effect of chi %%%
            v_chain2 = @(omega,omega_graft,chi) (v*sqrt((ones(size(omega))-2*chi*omega_graft).^2)).^(5/5);
            v_graft2 = v_chain2(kernel,zeros(size(omega_graft)),0);
            R_chain2 = @(omega,N,omega_graft,chi) ro*fo^(1/5)*(v_graft2).^(1/5).*beta_kuhn^(2/5)*N^(3/5).*omega.^(-3/5)*(beta_kuhn/ro)^(3/5) - 0.0*ro;
            E_chain2 = @(omega,N,omega_graft,chi) (-R_chain2(omega,N,omega_graft,chi).^2/(N*beta_kuhn^2) - v*fo*ro^2*N^2./((R_chain2(omega,N,omega_graft,chi)+ro).*omega).^3);            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            E_entropic_graft = E_entropic(kernel,N,omega_graft,chi);
            
            % Store grafting
            P_store{j} = P_graft;
            E_store{j} = E_graft;
            E_entropic_store{j} = E_entropic_graft;
            
            % Grabbing index
            p_indx = find(P_graft >= p_jth);
            
            % Select max
            if (j ~= 1) || (length(p_indx) > 1)
                P_indx_graft = P_graft(p_indx);
                [C_indx,C_ia] = setdiff(p_indx,graft_indx);
                [~,indx_graft] = max(P_graft(C_indx));
                graft_jth = C_ia(indx_graft);
                if isempty(graft_jth)
                    graft_jth = 1;
                end
            else
                % Pick random index from selected values
                graft_jth = randi([1,length(p_indx)]);
            end
            % Select open points
            if isempty(find(graft_indx == p_indx(graft_jth), 1))
                graft_indx(j) = p_indx(graft_jth);
                R_graft_tmp = R_chain(kernel,N,omega_graft,chi);
                flag_add = 1;
            end
        end
    end
    
    % Update R
    R_sum = zeros(size(R_graft_tmp));
    R_sum(graft_indx) = R_graft_tmp(graft_indx);
    R_graft = R_graft + R_sum;
    
    % Update graft counters
    for j = 1:n_graft
        graft_count(graft_indx(j)) = graft_count(graft_indx(j)) + 1;
    end
    
    % Store
    P_mc_store{ii} = P_store;
    E_mc_store{ii} = E_store;
    E_entropic_mc_store{ii} = E_entropic_store;
end

% Get distances
graft_count = graft_count/n_mc;
R_graft = R_graft/n_mc;

% Temp arrays
R_tmp = zeros(length(pts_poly(:,1)),1);
for j = 1:numel(R_tmp(:,1))
    R_tmp(j) = (R_graft(j)+ro - min(R_graft+ro));
end
% Scale distance
% R_tmp = (R_tmp - min(R_tmp(:)))/(max(R_tmp(:))-min(R_tmp(:)));
x_total = zeros(length(pts_poly(:,1)),1);
y_total = zeros(length(pts_poly(:,1)),1);
z_total = zeros(length(pts_poly(:,1)),1);
% Generate corona
for j = 1:numel(x_total)
    n_vect = [pts_poly(j,1) pts_poly(j,2) pts_poly(j,3)]/norm([pts_poly(j,1) pts_poly(j,2) pts_poly(j,3)]);
    x_total(j) = n_vect(1)*R_tmp(j) + pts_poly(j,1);
    y_total(j) = n_vect(2)*R_tmp(j) + pts_poly(j,2);
    z_total(j) = n_vect(3)*R_tmp(j) + pts_poly(j,3);
end

% Outputs
pts_core = pts_poly;
pts_corona = [x_total y_total z_total];

[theta,phi,r] = xyz_to_kernel(pts_core);

% Average all quantities
P_avg = cell(1,n_graft);
E_avg = cell(1,n_graft);
E_entropic_avg = cell(1,n_graft);
for nn = 1:n_graft
    P_tmp = zeros(size(P_store{1}));
    E_tmp = zeros(size(P_store{1}));
    E_entropic_tmp = zeros(size(P_store{1}));
    for i = 1:n_mc
        P_tmp = P_tmp + P_mc_store{i}{nn};
        E_tmp = E_tmp + E_mc_store{i}{nn};
        E_entropic_tmp = E_entropic_tmp + E_entropic_mc_store{i}{nn};
    end
    P_avg{nn} = P_tmp/n_mc;
    E_avg{nn} = E_tmp/n_mc;
    E_entropic_avg{nn} = E_entropic_tmp/n_mc;
end

theta_round = roundn(theta*180/pi,0);
theta_unique = unique(theta_round);

P_theta_avg = cell(1,n_graft);
E_theta_avg = cell(1,n_graft);
E_theta_entropic_avg = cell(1,n_graft);
for nn = 1:n_graft
    
    P_tmp = P_avg{nn};
    E_tmp = E_avg{nn};
    E_entropic_tmp = E_entropic_avg{nn};
    
    P_nn = zeros(size(theta_unique));
    E_nn = zeros(size(theta_unique));
    E_entropic_nn = zeros(size(theta_unique));
    for i = 1:length(theta_unique)
        indx_tmp = find(theta_round == theta_unique(i));
        P_nn(i) = mean(P_tmp(indx_tmp));
        E_nn(i) = mean(E_tmp(indx_tmp));
        E_entropic_nn(i) = mean(E_entropic_tmp(indx_tmp));
    end
    P_theta_avg{nn} = P_nn;
    E_theta_avg{nn} = E_nn;
    E_theta_entropic_avg{nn} = E_entropic_nn;
end

% Save data
mat_name = strcat(shape,'_sigma_',num2str(sigma,3),'_N_',num2str(N,3),'_kT_',num2str(kT,3),'.mat');
disp(mat_name)
save(mat_name)
