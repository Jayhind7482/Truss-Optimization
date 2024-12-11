%
% A program to generate the ground structure for truss topology optimization in 3D
%

% Clear workspace
clear all
clc

% User inputs
Lx = 8; % Length along x-axis
Ly = 4; % Length along y-axis
Lz = 2; % Length along z-axis
nNodex = 6; % Number of nodes along x-axis
nNodey = 3; % Number of nodes along y-axis
nNodez = 3; % Number of nodes along z-axis
E = 210E9; % Young's Modulus of the material

% Number of elements along each axis
Nx = nNodex - 1;
Ny = nNodey - 1;
Nz = nNodez - 1;

% Generate grid nodes
NNODE = 0;
for k = 1:Nz+1
    for j = 1:Ny+1
        for i = 1:Nx+1
            NNODE = NNODE + 1;
            nx(NNODE) = (Lx/Nx) * (i-1);
            ny(NNODE) = (Ly/Ny) * (j-1);
            nz(NNODE) = (Lz/Nz) * (k-1);
        end
    end
end

% Generate grid elements
NELEM = 0;
for i = 1:NNODE
    for j = i+1:NNODE
        NELEM = NELEM + 1;
        ncon(NELEM, 1) = i;
        ncon(NELEM, 2) = j;
    end
end

Ee = ones(NELEM, 1) * E;

% Compute element lengths
for i = 1:NELEM
    n1 = ncon(i, 1);
    n2 = ncon(i, 2);
    length(i) = sqrt((nx(n2)-nx(n1))^2 + (ny(n2)-ny(n1))^2 + (nz(n2)-nz(n1))^2);
end

%% Removing Overlapping Elements Based on Direction Cosines
tol = 1e-5; % Tolerance for identifying duplicate direction cosines
removed = true; % Flag to indicate if overlaps are removed in the current iteration

while removed
    removed = false; % Reset removal flag

    for node = 1:NNODE
        % Find all elements connected to the current node
        connected_elems = find(ncon(:, 1) == node | ncon(:, 2) == node);
        
        % Skip processing if no elements are connected
        if isempty(connected_elems)
            continue;
        end

        dir_cosines = []; % Store direction cosines of connected elements
        elem_to_keep = true(numel(connected_elems), 1); % Logical array for keeping elements

        % Compute direction cosines for all connected elements
        for i = 1:numel(connected_elems)
            elem = connected_elems(i);
            other_node = ncon(elem, 1) + ncon(elem, 2) - node; % Get the other node
            dx = nx(other_node) - nx(node);
            dy = ny(other_node) - ny(node);
            dz = nz(other_node) - nz(node);
            norm_val = sqrt(dx^2 + dy^2 + dz^2);
            dir_cosines(i, :) = [dx, dy, dz] / norm_val; % Direction cosines
        end

        % Compare direction cosines to find duplicates
        for i = 1:numel(connected_elems)
            for j = i+1:numel(connected_elems)
                if elem_to_keep(i) && elem_to_keep(j)
                    % Check if direction cosines are the same within tolerance
                    if all(abs(dir_cosines(i, :) - dir_cosines(j, :)) < tol)
                        % Remove the element with the larger length
                        if length(connected_elems(i)) > length(connected_elems(j))
                            elem_to_keep(i) = false;
                        else
                            elem_to_keep(j) = false;
                        end
                        removed = true; % Mark that an element was removed
                    end
                end
            end
        end

        % Remove unnecessary elements
        connected_elems_to_remove = connected_elems(~elem_to_keep);
        ncon(connected_elems_to_remove, :) = 0; % Mark as invalid
    end

    % Remove all marked elements
    valid_elems = all(ncon ~= 0, 2); % Find valid elements
    ncon = ncon(valid_elems, :);
    length = length(valid_elems);
    NELEM = size(ncon, 1);
end

%% Plotting the 3D Ground Structure
figure(1);
clf;
hold on;

% Plot elements
for i = 1:NELEM
    eye = ncon(i, 1);
    jay = ncon(i, 2);
    plot3([nx(eye), nx(jay)], [ny(eye), ny(jay)], [nz(eye), nz(jay)], '-r', 'LineWidth', 1.2);
end

% Plot nodes
for i = 1:NNODE
    plot3(nx(i), ny(i), nz(i), 'b*', 'MarkerSize', 8);
    text(nx(i), ny(i), nz(i), num2str(i), 'FontSize', 10, 'HorizontalAlignment', 'right');
end

% Configure plot
xlabel('X-axis', 'FontSize', 12);
ylabel('Y-axis', 'FontSize', 12);
zlabel('Z-axis', 'FontSize', 12);
title('3D Ground Structure for Truss Topology Optimization', 'FontSize', 16);
grid on;
axis equal;
view(3);
hold off;

% Save data
save trussdata3D NNODE NELEM nx ny nz ncon length Ee
