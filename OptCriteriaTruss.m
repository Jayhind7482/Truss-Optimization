clc;
clear all;
load('trussdata3D.mat');
tol = 1e-5;  % Tolerance for convergence
Le = length;
Ae = ones(1, NELEM);
Amax = 15; % Maximum allowed cross-sectional area of truss member
Amin = 1e-6; % Minimum allowed cross-sectional area of truss member
Ve = Ae .* Le;
vstar = 1000; % Maximum material available
Ae = Ae * vstar / sum(Le);

%% Force and Boundary Conditions (to be specified by the user)
%{

F(30, 1) = 20e4;
F(45, 1) = 00.20e4;
F(15, 1) = 0.2e4;
F(60, 1) = 2e4;
%}

%{
F(30, 1) = 0.5e4;
F(180, 1) = 0.5e4;
F(330, 1) = 0.5e4;
F(480, 1) = 0.5e4;
F(630, 1) = 0.5e4;
F(750, 1) = -0.5e4;
F(600, 1) = -0.5e4;
F(450, 1) = -0.5e4;
F(300, 1) = -0.5e4;
F(150, 1) = -0.5e4;

F(629, 1) = 0.5e4;
F(659, 1) = 0.5e4;
F(689, 1) = 0.5e4;
F(719, 1) = 0.5e4;
F(749, 1) = 0.5e4;
F(149, 1) = -0.5e4;
F(119, 1) = -0.5e4;
F(89, 1) = -0.5e4;
F(59, 1) = -0.5e4;
F(29, 1) = -0.5e4; 
%}


% F = zeros(3 * NNODE, 1);
%{
F(123, 1) = -50e3;
F(126, 1) = -50e3;
F(129, 1) = -50e3;
F(132, 1) = -50e3;
F(135, 1) = -50e3;
F(138, 1) = -50e3;
F(141, 1) = -50e3;%}
F(144, 1) = -50e3;
F(147, 1) = -50e3;
F(150, 1) = -50e3;
F(153, 1) = -50e3;
F(156, 1) = -50e3;
F(159, 1) = -50e3;
F(162, 1) = -50e3;
F(165, 1) = -50e3;
F(168, 1) = -50e3;
F(171, 1) = -50e3;
F(174, 1) = -50e3;
F(177, 1) = -50e3;
F(180, 1) = -50e3;
F(183, 1) = -50e3;
F(186, 1) = -50e3;
F(189, 1) = -50e3;
F(192, 1) = -50e3;
F(195, 1) = -50e3;
F(198, 1) = -50e3;
F(201, 1) = -50e3;
F(204, 1) = -50e3;
F(207, 1) = -50e3;
F(210, 1) = -50e3;
F(213, 1) = -50e3;
F(216, 1) = -50e3;
%}

F = zeros(3 * NNODE, 1);
F(114, 1) = 20e4;
F(117, 1) = 20e4;
F(120, 1) = 20e4;
F(123, 1) = 20e4;
F(132 ,1) = 20e4;
F(135 ,1) = 20e4;
F(138 ,1) = 20e4;
F(141 ,1) = 20e4;
F(150 ,1) = 20e4;
F(153 ,1) = 20e4;
F(156 ,1) = 20e4;
F(159 ,1) = 20e4;

%dispID = [1,2,3,16,17,18,31,32,33,43,44,45,241,242,243,256,257,258,271,272,273,286,287,288];
%dispID = [1, 2, 3, 31, 32, 33, 61, 62, 63, 91, 92, 93, 121, 122, ...
%        123, 151, 152, 153, 181, 182, 183, 211, 212, 213, 241, 242, 243, ... 
%       271, 272, 273, 301, 302, 303, 331, 332, 333, 361, 362, 363, 391, 392, .... 
%        393, 421, 422, 423, 451, 452, 453, 481, 482, 483, 511, 512, 513, 541, ...
%       542, 543, 571, 572, 573, 601, 602, 603, 631, 632, 633, 661, 662, 663, ... 
%        691, 692, 693, 721, 722, 723];

%{
dispID = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, ...
        20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, ...
        38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, ...
        56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,  ...
        74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, ...
        92, 93, 94, 95, 96];
%}


dispID = [1,2,3,16,17,18,19,20,21,34,35,36,37,38,39,52,53,54] ;

dispVal = zeros(numel(dispID), 1);

%% Optimality criteria method
beta = 0.008;    %(Tune)
its = 100; % Max number of iterations

% Initialize strain energy history
strain_energy_history = zeros(1, its);

history_A = zeros(NELEM, its + 1); % Store area values at each iteration
history_A(:, 1) = Ae(1, :);

for i = 1:its
    lambda1 = 0;
    [u, Rdisp, P, Ksing, SE, SEgrad] = femtruss3D(Ae, Le, Ee, nx, ny, nz, ncon, NELEM, NNODE, F, dispID);
    disp(SE);
    % Store total strain energy
    strain_energy_history(i) = SE +0.3; % Summing strain energy for all elements

    % Compute strain energy density
    for ie = 1:NELEM
        eye = ncon(ie, 1);
        jay = ncon(ie, 2);
        L = Le(ie);
        A = Ae(ie);
        V = A * L;
        E = Ee(ie);
        lx = (nx(jay) - nx(eye)) / L;
        ly = (ny(jay) - ny(eye)) / L;
        lz = (nz(jay) - nz(eye)) / L;

        % Transformation matrix
        Lambda = [lx, ly, lz,  0,  0,  0;
                   0,  0,  0, lx, ly, lz];

        % Local stiffness matrix
        k_bar = [1, -1; -1, 1];
        k1 = k_bar * (A * E / L);
        klocal = Lambda' * k1 * Lambda;

        % Element stiffness and strain energy
        Ke = klocal;
        ue = [u(3 * eye - 2); u(3 * eye - 1); u(3 * eye); ...
              u(3 * jay - 2); u(3 * jay - 1); u(3 * jay)];
        SEe = ue' * Ke * ue;
        SEVe(ie) = SEe / V;
    end

    % Compute new areas using optimality criteria
    for j = 1:NELEM
        lambda1 = lambda1 + Ae(j) * Le(j) * (SEVe(j))^beta;
    end
    lambda = (lambda1 / vstar)^(1 / beta);
    for j = 1:NELEM
        Ae(j) = Ae(j) * (SEVe(j) / lambda)^beta;
    end

    % Enforce area constraints
    flag = false;
    for j = 1:NELEM
        if Ae(j) > Amax
            Ae(j) = Amax;
            flag = true;
        end
        if Ae(j) < Amin
            Ae(j) = Amin;
            flag = true;
        end
    end

    while flag
        % Recompute areas under new constraints
        lambda1 = 0;
        Vup = 0;
        Vlow = 0;
        for j = 1:NELEM
            if Ae(j) == Amax
                Vup = Vup + Ae(j) * Le(j);
                continue;
            end
            if Ae(j) == Amin
                Vlow = Vlow + Ae(j) * Le(j);
                continue;
            end
            lambda1 = lambda1 + Ae(j) * Le(j) * (SEVe(j))^beta;
        end
        vstar1 = vstar - (Vup + Vlow);
        lambda = (lambda1 / vstar1)^(1 / beta);
        for j = 1:NELEM
            if Ae(j) ~= Amax && Ae(j) ~= Amin
                Ae(j) = Ae(j) * (SEVe(j) / lambda)^beta;
            end
        end
        flag = false;
        for j = 1:NELEM
            if Ae(j) > Amax
                Ae(j) = Amax;
                flag = true;
            end
            if Ae(j) < Amin
                Ae(j) = Amin;
                flag = true;
            end
        end
    end

    history_A(:, i + 1) = Ae(1, :);

    % Check for convergence
    tolcheck = max(abs(history_A(:, i + 1) - history_A(:, i)));
    if tolcheck <= tol
        break;
    end
end

Ae = (history_A(:, i + 1))';

% Filter and prepare optimized structure
j = 1;
for i = 1:NELEM
    if Ae(i) > 0.5
        NCON(j, :) = ncon(i, :);
        AE(j) = Ae(i);
        j = j + 1;
    end
end

%% PLOTING THE OPTIMIZED 3D TRUSS STRUCTURE
figure(1);
clf;
axis equal;
[NELEM1, NODE] = size(NCON);
for i = 1:NELEM1
    eye = NCON(i, 1);
    jay = NCON(i, 2);
    lw = AE(i);
    hold on;
    plot3([nx(eye), nx(jay)], [ny(eye), ny(jay)], [nz(eye), nz(jay)], '-k', 'LineWidth', lw);
end

% Formatting the plot
plot3(nx, ny, nz, 'k.', 'MarkerSize', 10); % Nodes
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Optimized 3D Truss Structure');
grid on;
view(3);

%% Plot Strain Energy vs Iteration
figure(2);
plot(1:its, strain_energy_history, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Iteration (i)', 'FontSize', 12);
ylabel('Total Strain Energy (SE)', 'FontSize', 12);
title('Strain Energy vs Iteration', 'FontSize', 14);
grid on;
