function [I_tensor,r_com,J_tensor,T0,v_I] = calc_inertial_tensor(pts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Assumes uniform density %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% Define faces %
%%%%%%%%%%%%%%%%
% Triangulate faces
k = boundary(pts,0);
% Sort by indices
[~,indx] = sort(sum(k.^2,2));
k = k(indx,:);

[~,mass] = boundary(pts,0);

% Define normal vector
n_store = cell(1,length(k(:,1)));
for i = 1:length(k(:,1))
    % Get face indices
    k_tmp = k(i,:);
    % Calculate normal vector
    % Displacement vector
    v1 = pts(k_tmp(2),:) - pts(k_tmp(1),:);
    v2 = pts(k_tmp(3),:) - pts(k_tmp(2),:);
    % Cross product
    v12 = cross(v1,v2);
    v12 = v12/norm(v12);
    % Storing normal vector
    n_store{i} = v12;
end

% Define face
face_store = cell(1,length(k(:,1)));
for i = 1:length(k(:,1))
    % Grab points
    k_tmp = [k(i,:) k(i,1)];
    pts_tmp = pts(k_tmp,:);
    % Define normal vector and normallize
    pts_norm = [mean(pts(k(i,:),:)); mean(pts(k(i,:),:)) + 0.5*n_store{i}];
    dn = pts_norm(2,:) - pts_norm(1,:);
    dn = dn/norm(dn);
    % Rotate to axis
    % Define axis
    z_new = dn;
    y_new = cross([1 0 0],z_new);
    x_new = cross(z_new,y_new);
    % Normalize
    y_new = y_new/norm(y_new);
    x_new = x_new/norm(x_new);
    % Rotate
    x_tmp = sum(bsxfun(@times,pts_tmp,x_new),2);
    y_tmp = sum(bsxfun(@times,pts_tmp,y_new),2);
    z_tmp = sum(bsxfun(@times,pts_tmp,z_new),2);
    % New rotated point
    pts_tmp = [x_tmp y_tmp z_tmp];
    % Determine angle
    angle_tmp = atan2(y_tmp,x_tmp);
    angle_tmp = angle_tmp(1:end-1);
    [~,indx] = sort(angle_tmp,'descend');
    % Set face order
    face_store{i} = k(i,indx) - 1;
end

% Update normal vector
n_store = cell(1,length(k(:,1)));
w_store = zeros(1,length(k(:,1)));
for i = 1:length(face_store)
    % Get face indices
    k_tmp = face_store{i} + 1;
    % Calculate normal vector
    % Displacement vector
    v1 = pts(k_tmp(2),:) - pts(k_tmp(1),:);
    v2 = pts(k_tmp(3),:) - pts(k_tmp(2),:);
    % Cross product
    v12 = cross(v1,v2);
    v12 = v12/norm(v12);
    % Offset from first point
    w12 = -1.0*dot(v12,pts(k_tmp(1),:));
    % Storing normal vector
    n_store{i} = v12;
    w_store(i) = w12;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate integrals for inertial tensor %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define working constants
X = 1;
Y = 2;
Z = 3;

% Initialize
T0 = 0;
T1 = zeros(1,3);
T2 = zeros(1,3);
TP = zeros(1,3);

% Define differnt masses to face
m_face = ones(1,length(face_store));

% Loop through faces
for i = 1:length(face_store)
    %%%%%%%%%%%%%%%
    % Select face %
    face_tmp = face_store{i} + 1;
    % Select face normal
    n_tmp = n_store{i};
    % Mass of polyhedron described by face
    m_i = m_face(i);
    %%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%
    % Define constants
    nx = abs(n_tmp(1));
    ny = abs(n_tmp(2));
    nz = abs(n_tmp(3));
    % Define x,y,z indices based on normal vector
    if (nx > ny) && (nx > nz)
        C = 0;
    else
        if (ny > nz)
            C = 1;
        else
            C = 2;
        end
    end
    A = mod(C+1,3);
    B = mod(A+1,3);
    % Update (index counting)
    % A B C
    A = A + 1;
    B = B + 1;
    C = C + 1;
    %%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate projection integrals
    % Initialize
    P1 = 0;
    Pa = 0;
    Pb = 0;
    Paa = 0;
    Pab = 0;
    Pbb = 0;
    Paaa = 0;
    Paab = 0;
    Pabb = 0;
    Pbbb = 0;
    for j = 1:length(face_tmp)
        % First set of vector
        a0 = pts(face_tmp(j),A);
        b0 = pts(face_tmp(j),B);
        % Second set of vector
        indx_tmp = j + 1;
        if j == length(face_tmp)
            indx_tmp = 1;
        end
        a1 = pts(face_tmp(indx_tmp),A);
        b1 = pts(face_tmp(indx_tmp),B);
        % Vector difference
        da = a1 - a0;
        db = b1 - b0;
        % Working parameters
        a0_2 = a0*a0;
        a0_3 = a0_2*a0;
        a0_4 = a0_3*a0;
        b0_2 = b0*b0;
        b0_3 = b0_2*b0;
        b0_4 = b0_3*b0;
        a1_2 = a1*a1;
        a1_3 = a1_2*a1;
        b1_2 = b1*b1;
        b1_3 = b1_2*b1;
        
        % Projection parameters
        C1 = a1 + a0;
        Ca = a1*C1 + a0_2;
        Caa = a1*Ca + a0_3;
        Caaa = a1*Caa + a0_4;
        Cb = b1*(b1 + b0) + b0_2;
        Cbb = b1*Cb + b0_3;
        Cbbb = b1*Cbb + b0_4;
        Cab = 3*a1_2 + 2*a1*a0 + a0_2;
        Kab = a1_2 + 2*a1*a0 + 3*a0_2;
        Caab = a0*Cab + 4*a1_3;
        Kaab = a1*Kab + 4*a0_3;
        Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
        Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;
        
        % Projection integrals
        P1 = P1 + db*C1;
        Pa = Pa + db*Ca;
        Paa = Paa + db*Caa;
        Paaa = Paaa + db*Caaa;
        Pb = Pb + da*Cb;
        Pbb = Pbb + da*Cbb;
        Pbbb = Pbbb + da*Cbbb;
        Pab = Pab + db*(b1*Cab + b0*Kab);
        Paab = Paab + db*(b1*Caab + b0*Kaab);
        Pabb = Pabb + da*(a1*Cabb + a0*Kabb);
    end
    % Prefactors
    P1 = P1/2;
    Pa = Pa/6;
    Paa = Paa/12;
    Paaa = Paaa/20;
    Pb = -1.0*Pb/6;
    Pbb = -1.0*Pbb/12;
    Pbbb = -1.0*Pbbb/20;
    Pab = Pab/24;
    Paab = Paab/60;
    Pabb = -1.0*Pabb/60;
    
    % Face integral constants
    w = w_store(i);
    k1 = 1.0/n_tmp(C); 
    k2 = k1*k1; 
    k3 = k2*k1; 
    k4 = k3*k1;
    
    % Compute face integrals
    Fa = k1*Pa;
    Fb = k1*Pb;
    Fc = -k2*(n_tmp(A)*Pa + n_tmp(B)*Pb + w*P1);
    
    Faa = k1*Paa;
    Fbb = k1*Pbb;
    Fcc = k3*(n_tmp(A)*n_tmp(A)*Paa + 2.0*n_tmp(A)*n_tmp(B)*Pab + ...
        n_tmp(B)*n_tmp(B)*Pbb + w*(2.0*(n_tmp(A)*Pa + n_tmp(B)*Pb) + w*P1));
    
    Faaa = k1*Paaa;
    Fbbb = k1*Pbbb;
    Fccc = -k4*(n_tmp(A)*n_tmp(A)*n_tmp(A)*Paaa + ...
        3.0*n_tmp(A)*n_tmp(A)*n_tmp(B)*Paab + 3.0*n_tmp(A)*n_tmp(B)*n_tmp(B)*Pabb ...
        + n_tmp(B)*n_tmp(B)*n_tmp(B)*Pbbb + 3.0*w*(n_tmp(A)*n_tmp(A)*Paa + ...
        2.0*n_tmp(A)*n_tmp(B)*Pab + n_tmp(B)*n_tmp(B)*Pbb) + w*w*(3.0*(n_tmp(A)*Pa ...
        + n_tmp(B)*Pb) + w*P1));
    
    Faab = k1*Paab;
    Fbbc = -k2*(n_tmp(A)*Pabb + n_tmp(B)*Pbbb + w*Pbb);
    Fcca = k3*(n_tmp(A)*n_tmp(A)*Paaa + 2.0*n_tmp(A)*n_tmp(B)*Paab + ...
        n_tmp(B)*n_tmp(B)*Pabb + w*(2.0*(n_tmp(A)*Paa + n_tmp(B)*Pab) + w*Pa));
    
    % Calculate mass volume integral
    % T0
    if (A == X)
        T0 = T0 + n_tmp(X)*Fa;
    else
        if (B == X)
            T0 = T0 + n_tmp(X)*Fb;
        else
            T0 = T0 + n_tmp(X)*Fc;
        end
    end
    
    % Calculatge other volume integrals
    % T1
    T1(A) = T1(A) + m_i*n_tmp(A)*Faa;
    T1(B) = T1(B) + m_i*n_tmp(B)*Fbb;
    T1(C) = T1(C) + m_i*n_tmp(C)*Fcc;
    % T2
    T2(A) = T2(A) + m_i*n_tmp(A)*Faaa;
    T2(B) = T2(B) + m_i*n_tmp(B)*Fbbb;
    T2(C) = T2(C) + m_i*n_tmp(C)*Fccc;
    % TP
    TP(A) = TP(A) + m_i*n_tmp(A)*Faab;
    TP(B) = TP(B) + m_i*n_tmp(B)*Fbbc;
    TP(C) = TP(C) + m_i*n_tmp(C)*Fcca;
end

% Final integrals
T1 = T1/2;
T2 = T2/3;
TP = TP/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define inertial tensor %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% T0 = mass;

% Define center of mass
r = T1/T0;
% Inertial tensor
% Diagonal
J(X,X) = (T2(Y) + T2(Z));
J(Y,Y) = (T2(X) + T2(Z));
J(Z,Z) = (T2(X) + T2(Y));
% Cross terms
% XY
J(X,Y) = -TP(X);
J(Y,X) = -TP(X);
% XZ
J(X,Z) = -TP(Z);
J(Z,X) = -TP(Z);
% YZ
J(Y,Z) = -TP(Y);
J(Z,Y) = -TP(Y);

% Shift to center of mass
% Diagonal
J(X,X) = J(X,X) - T0*(r(Y)*r(Y) + r(Z)*r(Z));
J(Y,Y) = J(Y,Y) - T0*(r(X)*r(X) + r(Z)*r(Z));
J(Z,Z) = J(Z,Z) - T0*(r(Y)*r(Y) + r(X)*r(X));
% Cross terms
% XY
J(X,Y) = J(X,Y) + T0*r(X)*r(Y);
J(Y,X) = J(Y,X) + T0*r(X)*r(Y);
% XZ
J(X,Z) = J(X,Z) + T0*r(X)*r(Z);
J(Z,X) = J(Z,X) + T0*r(X)*r(Z);
% YZ
J(Y,Z) = J(Y,Z) + T0*r(Y)*r(Z);
J(Z,Y) = J(Z,Y) + T0*r(Y)*r(Z);

% Diagonalize to determine principle axis
[v_I,I_p] = eig(J);
I_p = I_p/(T0);

% Ouput variables
r_com = r;
I_tensor = I_p;
J_tensor = J/T0;
