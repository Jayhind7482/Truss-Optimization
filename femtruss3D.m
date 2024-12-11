function [u, Rdisp, P, Ksing, SE, SEgrad] = ...
         femtruss3D(Ae, Le, Ee, nx, ny, nz, ncon, NELEM, NNODE, F, dispID)

K = zeros(3 * NNODE, 3 * NNODE); % Initialize global stiffness matrix
SEgrad = zeros(NELEM, 1);

for ie = 1:NELEM
    eye = ncon(ie, 1);
    jay = ncon(ie, 2);

    % Compute transformation matrix, Lambda
    L = Le(ie);
    A = Ae(ie);
    E = Ee(ie);
    lx = (nx(jay) - nx(eye)) / L;
    ly = (ny(jay) - ny(eye)) / L;
    lz = (nz(jay) - nz(eye)) / L;

    Lambda = [lx, ly, lz,  0,  0,  0;
               0,  0,  0, lx, ly, lz];

    % Local stiffness matrix in local coordinates
    k = [1, -1; -1, 1] * (A * E / L);

    % Transform to global coordinates
    klocal = Lambda' * k * Lambda;

    % Assemble into global stiffness matrix
    id1 = 3 * (eye - 1) + 1;
    id2 = id1 + 1;
    id3 = id1 + 2;
    id4 = 3 * (jay - 1) + 1;
    id5 = id4 + 1;
    id6 = id4 + 2;

    ids = [id1, id2, id3, id4, id5, id6];

    for i = 1:6
        for j = 1:6
            K(ids(i), ids(j)) = K(ids(i), ids(j)) + klocal(i, j);
        end
    end
end

% Store K before applying boundary conditions
Ksing = K;

% Apply displacement boundary conditions
[sm, sn] = size(dispID);
Ndbc = sn;
for nd = 1:Ndbc
   
    K = matcut(K, dispID(nd) - nd + 1); % Remove rows/columns in K
    F = veccut(F, dispID(nd) - nd + 1);                % Remove corresponding entries in F
end

% Debug dimensions before solving
if size(K, 1) ~= length(F)
    error('Mismatch in dimensions of K (%d x %d) and F (%d)', size(K, 1), size(K, 2), length(F));
end

% Solve for unknown displacements
F;
U = inv(K)*F;

SE = 0.5 * U' * K * U;

% Store displacements including boundary conditions
u = zeros(3 * NNODE, 1);

for iu = 1:Ndbc
    u(dispID(iu)) = 12345.12345;
end

iuc = 0;
for iu = 1:3 * NNODE
    if u(iu) == 12345.12345
        iuc = iuc + 1;
    else
        u(iu) = U(iu - iuc);
    end
end

for iu = 1:Ndbc
    u(dispID(iu)) = 0;
end

dx = zeros(1, NNODE);
dy = zeros(1, NNODE);
dz = zeros(1, NNODE);

for iu = 1:NNODE
    dx(iu) = u(3 * (iu - 1) + 1);
    dy(iu) = u(3 * (iu - 1) + 2);
    dz(iu) = u(3 * (iu - 1) + 3);
end

% Compute reaction forces
R = Ksing * u;
Rdisp = zeros(1, Ndbc);
for iu = 1:Ndbc
    Rdisp(iu) = R(dispID(iu));
end

% Compute internal forces and strains
P = zeros(NNODE, 6);
for ie = 1:NELEM
    eye = ncon(ie, 1);
    jay = ncon(ie, 2);

    L = Le(ie);
    A = Ae(ie);
    lx = (nx(jay) - nx(eye)) / L;
    ly = (ny(jay) - ny(eye)) / L;
    lz = (nz(jay) - nz(eye)) / L;

    Lambda = [lx, ly, lz,  0,  0,  0;
               0,  0,  0, lx, ly, lz];

    k = [1, -1; -1, 1] * (A * E / L);

    klocal = Lambda' * k * Lambda;

    id1 = 3 * (eye - 1) + 1;
    id2 = id1 + 1;
    id3 = id1 + 2;
    id4 = 3 * (jay - 1) + 1;
    id5 = id4 + 1;
    id6 = id4 + 2;

    ulocal = [u(id1), u(id2), u(id3), u(id4), u(id5), u(id6)];

    Rint = klocal * ulocal';
    SEgrad(ie) = -0.5 * ulocal * klocal * ulocal' / A;

    P(ie, :) = Rint;
end
