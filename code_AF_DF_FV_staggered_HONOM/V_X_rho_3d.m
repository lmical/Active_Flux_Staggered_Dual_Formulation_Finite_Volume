clc;
close all;
clear all;

% Folder and filenames
folder = 'simulations_paper/sod2d_400/space_reconstruction-27/time_scheme_3/PP12/CFL0.475/';
mesh_x_file = fullfile(folder, 'W_X_MESH_CoordinateX.dat');
mesh_y_file = fullfile(folder, 'W_X_MESH_CoordinateY.dat');
solution_file = fullfile(folder, 'W_X_SOLUTION_OCT_0000000.2500000.dat');

% Parameters
delimiterIn = ' ';
headerlinesIn = 1;

% Check that files exist
assert(isfile(mesh_x_file), ['File not found: ', mesh_x_file]);
assert(isfile(mesh_y_file), ['File not found: ', mesh_y_file]);
assert(isfile(solution_file), ['File not found: ', solution_file]);

% Load mesh data
x = importdata(mesh_x_file, delimiterIn, headerlinesIn).data(:,1);
y = importdata(mesh_y_file, delimiterIn, headerlinesIn).data(:,1);
[X, Y] = meshgrid(x, y);

% Store domain limits
xIni = min(x); xEnd = max(x);
yIni = min(y); yEnd = max(y);

% Load solution data
sol = importdata(solution_file, delimiterIn, headerlinesIn).data;
ro = sol(:,1);
u  = sol(:,2);
v  = sol(:,3);
p  = sol(:,4);

% Reshape variables
[m, n] = size(X);
RO = reshape(ro, [n, m])';
U  = reshape(u , [n, m])';
V  = reshape(v , [n, m])';
P  = reshape(p , [n, m])';

% Compute derived quantities
Gmm = 1.4;
ROU = RO .* U;
ROV = RO .* V;
ENERGY = P / (Gmm - 1) + 0.5 * RO .* (U.^2 + V.^2);

% Plotting
figure('Position', [100, 100, 800, 400])
surf(X, Y, RO, RO, 'EdgeColor', 'none');
colormap('parula');
colorbar;
title('\rho', 'FontWeight', 'normal');
%xlabel('x');
%ylabel('y');
view(30, 30) % View angle, adjust as needed
set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1])
set(gca, 'YTick', [-1, -0.5, 0, 0.5, 1])
xlim([xIni, xEnd]);
ylim([yIni, yEnd]);
axis tight;
daspect([4 4 1.5]) % Adjust aspect ratio for clarity

cb = colorbar;
% Shrink the axes and move the colorbar closer
ax = gca;
ax_pos = ax.Position;
ax.Position = [ax_pos(1), ax_pos(2), ax_pos(3)*0.92, ax_pos(4)];  % Less shrinking
cb.Position = [ax.Position(1) + ax.Position(3) - 0.125, ax.Position(2), 0.02, ax.Position(4)];


% Save figure with high quality (HQ)
saveas(gcf, 'sod2d_V_X_3d_density.pdf') % Save as vector PDF (HQ)
% Or save as PNG with high resolution (if preferred)
%print(gcf, 'sod2d_U_3d_density.png', '-dpng', '-r300') % 300 DPI for high-res PNG
