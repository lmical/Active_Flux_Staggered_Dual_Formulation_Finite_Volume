clc
close all
clear all


% Parameters
folder = 'simulations_paper/sod2d_400/space_reconstruction-27/time_scheme_3/PP12/CFL0.475/';
filename = 'SOLUTION_OCT_0000000.2500000.dat';

% File paths
mesh_x_file = fullfile(folder, 'MESH_CoordinateX.dat');
mesh_y_file = fullfile(folder, 'MESH_CoordinateY.dat');
solution_file = fullfile(folder, filename);

% Check file existence
assert(isfile(mesh_x_file), ['File not found: ', mesh_x_file]);
assert(isfile(mesh_y_file), ['File not found: ', mesh_y_file]);
assert(isfile(solution_file), ['File not found: ', solution_file]);

% Load mesh data (skip header line)
x = dlmread(mesh_x_file, ' ', 1, 0);
y = dlmread(mesh_y_file, ' ', 1, 0);
[X, Y] = meshgrid(x, y);

% Load solution data (skip header line)
solution_data = dlmread(solution_file, '', 1, 0);
ro = solution_data(:, 1);
rou = solution_data(:, 2);
rov = solution_data(:, 3);
energy = solution_data(:, 4);

% Reshape and compute derived quantities
Gmm = 1.4;
[m, n] = size(X);
RO = reshape(ro, m, n);
ROU = reshape(rou, m, n);
ROV = reshape(rov, m, n);
ENERGY = reshape(energy, m, n);

% Avoid division by zero
U = zeros(size(RO));
V = zeros(size(RO));
U(RO ~= 0) = ROU(RO ~= 0) ./ RO(RO ~= 0);
V(RO ~= 0) = ROV(RO ~= 0) ./ RO(RO ~= 0);
P = (Gmm - 1) * (ENERGY - 0.5 * RO .* (U.^2 + V.^2));

% 3D surface plot
figure('Position', [100, 100, 800, 400])
surf(X, Y, RO, RO, 'EdgeColor', 'none')
colormap('parula') % You can change this to another colormap if needed
colorbar
%xlabel('x')
%ylabel('y')
%zlabel('\rho')
title('\rho', 'FontWeight', 'normal');
view(30, 30) % View angle, adjust as needed
set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1])
set(gca, 'YTick', [-1, -0.5, 0, 0.5, 1])
axis tight
daspect([4 4 1.5]) % Adjust aspect ratio for clarity

cb = colorbar;
% Shrink the axes and move the colorbar closer
ax = gca;
ax_pos = ax.Position;
ax.Position = [ax_pos(1), ax_pos(2), ax_pos(3)*0.92, ax_pos(4)];  % Less shrinking
cb.Position = [ax.Position(1) + ax.Position(3) - 0.125, ax.Position(2), 0.02, ax.Position(4)];

% Save figure with high quality (HQ)
saveas(gcf, 'sod2d_U_3d_density.pdf') % Save as vector PDF (HQ)
% Or save as PNG with high resolution (if preferred)
%print(gcf, 'sod2d_U_3d_density.png', '-dpng', '-r300') % 300 DPI for high-res PNG
