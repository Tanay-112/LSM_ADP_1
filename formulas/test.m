clear; clc;

params.b = 20;        % [mm]
params.c = 16;        % [mm]
params.f = 2;         % [mm]
params.d = params.c + params.f;
params.Ef = 70000;    % [N/mm^2]
params.E0 = 70000;    % [N/mm^2]
params.K = 0.5;
params.L = 100;       % [mm]
params.d_strut = 0.4; % [mm]
params.a_cell = 4;    % [mm]

[Pcr_global, mode_global, details_global] = global_formulas(params);
[Pcr_wrinkling, details_wrinkling] = wrinkling_formulas(params);

fprintf('Global buckling:\n');
fprintf('  Pcr_global  = %.4e N\n', Pcr_global);
fprintf('  mode_global = %s\n', mode_global);

fprintf('\n');
fprintf('Wrinkling:\n');
fprintf('  Pcr_wrinkling = %.4e N\n', Pcr_wrinkling);
fprintf('\n');