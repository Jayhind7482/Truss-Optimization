function [obj, grad] = truss_obj_grad(Area)
global Le F dispID Ee nx ny nz ncon NELEM NNODE Ve

% Compute displacements and strain energy
[u, Rdisp, P, Ksing, SE, SEgrad] = femtruss3D(Area, Le, Ee, nx, ny, nz, ncon, NELEM, NNODE, F, dispID);

% Objective function: compliance (F' * u)
obj = F' * u;

% Gradient of the objective function
grad = zeros(1, NELEM);
for ie = 1:NELEM
    eye = ncon(ie, 1); % Node i
    jay = ncon(ie, 2); % Node j

    % Element geometry and material properties
    L = Le(ie);
    E = Ee(ie);

    % Direction cosines
    lx = (nx(jay) - nx(eye)) / L;
    ly = (ny(jay) - ny(eye)) / L;
    lz = (nz(jay) - nz(eye)) / L;

    % Transformation matrix (3D)
    Lambda = [lx, ly, lz,  0,  0,  0;
               0,  0,  0, lx, ly, lz];

    % Local stiffness matrix
    k_bar = [1, -1; -1, 1];
    k = k_bar * (E / L);

    % Transform to global coordinates
    ke = Lambda' * k * Lambda;

    % Displacement vector for the element
    ue = [u(3 * eye - 2); u(3 * eye - 1); u(3 * eye); ...
          u(3 * jay - 2); u(3 * jay - 1); u(3 * jay)];

    % Gradient of the objective function
    grad(1, ie) = -ue' * ke * ue;
end
