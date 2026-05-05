%% ============================================================
%  SANDWICH PANEL BUCKLING ANALYSIS
%  Reproduces Fig. 10 from the reference paper:
%  "Buckling of lattice core sandwich panels under uniaxial
%   compressive load" (2D BCC lattice core, simply supported)
%
%  Two failure modes are analyzed:
%    1. Global buckling  -> FSDT-based analytical model (Eq. 23-27)
%    2. Local strut buckling -> Euler buckling with end constraint
%                              factor n (Eq. 28-31)
%
%  Critical load = min(F_local, F_global) at each slenderness ratio
%% ============================================================
clear; clc; close all;

%% ---------------------------------------------------------------
%  SECTION 1: MATERIAL PROPERTIES
%  Both face sheets and struts are assumed to be the same
%  aluminum alloy: Es = Ef = 70 GPa, nu = 0.35
% ----------------------------------------------------------------
Es   = 70000;           % Strut (solid) elastic modulus [MPa]
Ef   = 70000;           % Face sheet elastic modulus [MPa]
nu_s = 0.35;            % Poisson's ratio (both materials)
Gf   = Ef/(2*(1+nu_s)); % Face sheet shear modulus [MPa]

%% ---------------------------------------------------------------
%  SECTION 2: LATTICE CORE GEOMETRY
%  The core is a 2D BCC lattice with Nz layers through thickness.
%  Unit cell size: a = hc/Nz
%  Stiffness ratio: Ef/Exx_c = 1000 (parametric study in paper)
% ----------------------------------------------------------------
hc = 21;                % Core thickness [mm]
Nz = 7;                 % Number of unit cells through core thickness
a  = hc/Nz;             % Unit cell size [mm] -> a = 3 mm

stiffness_ratio = 1000; % Ef / Exx_c (parametric input)
Exx_c = Ef/stiffness_ratio; % Effective core modulus in load direction [MPa]

% --- Strut diameter from Eq. (1): Exx_c/Es = (sqrt(2)-1)*pi/4 * (d/a)^2 ---
C_Exx = (sqrt(2)-1)*pi/4;          % Analytical coefficient for Exx_c
d     = a*sqrt(Exx_c/(Es*C_Exx));  % Strut diameter [mm]

% --- Other effective core properties from Eq. (1) ---
Ezz_c = Es*(pi/4)          *(d/a)^2; % Effective modulus in thickness dir. [MPa]
Gxz_c = Es*(pi/(4*sqrt(2)))*(d/a)^2; % Effective shear modulus [MPa]

% --- Strut cross-section and lengths ---
Is   = pi*d^4/64;    % Second moment of area (circular cross-section) [mm^4]
ls_i = a/sqrt(2);    % Inclined strut length [mm] (45 deg in unit cell)
ls_v = a;            % Vertical strut length [mm]

fprintf('============================================================\n')
fprintf(' MATERIAL & GEOMETRY\n')
fprintf('============================================================\n')
fprintf('Unit cell size     : a     = %.3f mm\n', a)
fprintf('Strut diameter     : d     = %.5f mm\n', d)
fprintf('Second mom. inertia: Is    = %.4e mm^4\n', Is)
fprintf('Inclined strut len : ls_i  = %.4f mm\n', ls_i)
fprintf('Vertical strut len : ls_v  = %.4f mm\n', ls_v)
fprintf('Eff. core modulus  : Exx_c = %.4f MPa\n', Exx_c)
fprintf('Eff. core modulus  : Ezz_c = %.4f MPa\n', Ezz_c)
fprintf('Eff. shear modulus : Gxz_c = %.4f MPa\n\n', Gxz_c)

%% ---------------------------------------------------------------
%  SECTION 3: END CONSTRAINT FACTOR n
%
%  The strut ends are neither pinned nor clamped -- they are
%  elastically restrained by neighboring struts (rotational spring).
%  Two methods are used to find n:
%
%  Method 1 (Eq. 38): Determinant of boundary condition matrix = 0
%    -> Solve: (3 + lambda^2*ls_i*ls_v)*sin(ls_i*lambda)
%              - 3*lambda*ls_i*cos(ls_i*lambda) = 0
%    -> n = ls_i*lambda / pi
%
%  Method 2 (Eq. 47): Stability stiffness matrix equilibrium
%    -> Solve: -2*s*c^2 + 2*s + s1/sqrt(2) + s1*c1/sqrt(2) = 0
%    where s,c depend on inclined strut (compression)
%          s1,c1 depend on vertical strut (tension)
% ----------------------------------------------------------------

%--- Method 1: Determinant condition (Eq. 38) ---
% Substitution x = ls_i * lambda to find the smallest positive root
eq38 = @(x)(3 + (x/ls_i)^2*ls_i*ls_v)*sin(x) - 3*(x/ls_i)*ls_i*cos(x);
x_scan = linspace(0.5, 8, 5000);
f_scan = arrayfun(eq38, x_scan);
idx1   = find(diff(sign(f_scan)) ~= 0, 1, 'first'); % find sign change
n_det  = fzero(eq38, [x_scan(idx1), x_scan(idx1+1)]) / pi;

%--- Method 2: Stiffness matrix equilibrium (Eq. 47) ---
eq47  = @(n) stiffness_eq47(n, ls_i, ls_v, Es, Is);
n_scan = linspace(0.6, 2.0, 2000);
f_scan2 = arrayfun(eq47, n_scan);
idx2    = find(diff(sign(f_scan2)) ~= 0, 1, 'first');
n_stiff = fzero(eq47, [n_scan(idx2), n_scan(idx2+1)]);

fprintf('============================================================\n')
fprintf(' END CONSTRAINT FACTOR n\n')
fprintf('============================================================\n')
fprintf('n (determinant method,    Eq.38) = %.4f  [paper: ~1.15]\n', n_det)
fprintf('n (stiffness matrix meth, Eq.47) = %.4f  [paper: ~1.17]\n\n', n_stiff)

%% ---------------------------------------------------------------
%  SECTION 4: LOCAL STRUT BUCKLING LOAD (Eq. 28-31)
%
%  Step 1 - Euler buckling load of one inclined strut (Eq. 31):
%    F_buck = n^2 * pi^2 * Es * Is / ls_i^2   [N]
%
%  Step 2 - Unit cell load at buckling:
%    The unit cell contains 4 inclined struts.
%    Under uniaxial load F [N/mm], the paper defines:
%      F_cell,i = -F_cell / sqrt(2)  (load on inclined struts, Eq. 30)
%    where F_cell,i represents 2 struts in compression.
%    -> Each single strut carries: F_cell / (2*sqrt(2))
%    -> At buckling: F_cell / (2*sqrt(2)) = F_buck
%    -> F_cell_crit = 2*sqrt(2) * F_buck
%
%  Step 3 - Convert unit cell load back to sandwich load F (Eq. 30 inverted):
%    From Eq. (29)-(30): F_cell = F * a * Exx_c / (Exx_c*hc + 2*Ef*hf)
%    Inverting:
%      F_local = F_cell_crit / a * (Exx_c*hc + 2*Ef*hf) / (Exx_c*hc)
%
%  Note: F_local is independent of sandwich length l (constant plateau)
% ----------------------------------------------------------------

% Euler buckling load for one inclined strut [N]
F_buck_det   = n_det^2   * pi^2 * Es * Is / ls_i^2;
F_buck_stiff = n_stiff^2 * pi^2 * Es * Is / ls_i^2;

thickness_ratios = [14, 21, 28]; % hc/hf values to analyze
ntr = length(thickness_ratios);

F_local_det   = zeros(1, ntr);
F_local_stiff = zeros(1, ntr);

fprintf('============================================================\n')
fprintf(' LOCAL STRUT BUCKLING LOADS\n')
fprintf('============================================================\n')
for i = 1:ntr
    hf = hc / thickness_ratios(i);
    axial_stiff = Exx_c*hc + 2*Ef*hf;  % Total axial stiffness [N/mm]

    % Critical sandwich load [N/mm] -- from inverted Eq. (30)
    F_local_det(i)   = 2*sqrt(2)*F_buck_det   / a * axial_stiff/(Exx_c*hc);
    F_local_stiff(i) = 2*sqrt(2)*F_buck_stiff / a * axial_stiff/(Exx_c*hc);

    fprintf('hc/hf = %2d | hf = %.3f mm | F_local(n=%.4f) = %7.2f N/mm',...
        thickness_ratios(i), hf, n_det, F_local_det(i))
    fprintf(' | F_local(n=%.4f) = %7.2f N/mm\n',...
        n_stiff, F_local_stiff(i))
end
fprintf('\n')

%% ---------------------------------------------------------------
%  SECTION 5: GLOBAL BUCKLING LOAD (FSDT, Eq. 23-27)
%
%  First-order Shear Deformation Theory (FSDT) is used.
%  The sandwich is modeled as a shear-deformable beam.
%
%  Key quantities:
%    k  = shear correction factor (Eq. 9-11)
%    S  = effective shear stiffness = k*(2*Gf*hf + Gxz_c*hc) [N/mm]
%    D  = bending stiffness [N.mm]
%
%  Critical load for mode m (Eq. 27):
%    F_cr(m) = (T12^2 + t11*T22) / (m*pi/l)^2 / T22)
%
%  Minimum over m = 1,2,3,... gives the critical global load.
% ----------------------------------------------------------------

slenderness = linspace(2, 60, 400); % l/hc range
colors_hex  = {'#1f77b4', '#ff7f0e', '#2ca02c'};

figure('Position', [80 80 720 540], 'Color', 'w');
hold on; box on; grid on;

fprintf('============================================================\n')
fprintf(' GLOBAL BUCKLING -- SHEAR STIFFNESS AND BENDING STIFFNESS\n')
fprintf('============================================================\n')

for i = 1:ntr
    hf = hc/thickness_ratios(i);

    % --- Shear correction factor k (Eq. 9-11) ---
    % Voigt (upper bound) and Reuss (lower bound) averages of shear modulus
    G_Voigt    = (2*Gf*hf + Gxz_c*hc) / (2*hf + hc);
    inv_G_Reuss = 2*hf/(Gf*(2*hf+hc)) + hc/(Gxz_c*(2*hf+hc));
    k           = (1/inv_G_Reuss) / G_Voigt;

    % --- Effective sandwich stiffnesses ---
    S = k*(2*Gf*hf + Gxz_c*hc);   % Shear stiffness [N/mm]
    D = (1/12)*(Exx_c*hc^3 + 2*Ef*hf*(4*hf^2 + 6*hf*hc + 3*hc^2)); % Bending [N.mm]

    fprintf('hc/hf = %2d | k = %.5f | S = %8.3f N/mm | D = %.4e N.mm\n',...
        thickness_ratios(i), k, S, D)

    % --- Compute global critical load for each slenderness ratio ---
    F_global = zeros(size(slenderness));
    for j = 1:length(slenderness)
        l = slenderness(j)*hc;   % Sandwich length [mm]
        Fcr_min = inf;

        % Sweep half-wave numbers m = 1..15 and take minimum
        for m = 1:15
            km  = m*pi/l;
            t11 = S*km^2;            % (Eq. 24)
            T12 = -S*km;             % (Eq. 25)
            T22 = -S - D*km^2;       % (Eq. 26)
            Fcr = (T12^2 + t11*T22) / (km^2 * T22); % (Eq. 27)
            if Fcr > 0 && Fcr < Fcr_min
                Fcr_min = Fcr;
            end
        end
        F_global(j) = Fcr_min;
    end

    % --- Critical load = min(local, global) at each slenderness ---
    F_cr_det   = min(F_local_det(i)   * ones(size(slenderness)), F_global);
    F_cr_stiff = min(F_local_stiff(i) * ones(size(slenderness)), F_global);

    % --- Plot ---
    rgb = hex2rgb(colors_hex{i});
    plot(slenderness, F_cr_stiff, '-',  'Color', rgb, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('h^{(c)}/h^{(f)} = %d', thickness_ratios(i)));
    plot(slenderness, F_cr_det,   '--', 'Color', rgb, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
end

% Add legend entries for line styles
plot(nan, nan, 'k-',  'LineWidth', 2,   'DisplayName', sprintf('Present model n = %.2f', n_stiff))
plot(nan, nan, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('Present model n = %.2f', n_det))

%% ---------------------------------------------------------------
%  SECTION 6: FIGURE FORMATTING
% ----------------------------------------------------------------
xlabel('Slenderness ratio $l/h^{(c)}$',   'Interpreter', 'latex', 'FontSize', 13)
ylabel('Critical load $F_{cr}$ in N/mm',  'Interpreter', 'latex', 'FontSize', 13)
title({'Sandwich Panel Buckling -- Fig.10 Reproduce', ...
    '$E^{(f)}/E^{(c)}_{xx} = 1000$,  $N_z = 7$,  $a = 3$ mm'}, ...
    'Interpreter', 'latex', 'FontSize', 11)
legend('Location', 'northeast', 'FontSize', 9)
xlim([0 60]);  ylim([0 1400])
xticks(0:10:60);  yticks(0:200:1400)
set(gca, 'FontSize', 11)
xline(17, ':', 'Color', [.5 .5 .5], 'LineWidth', 1.2, ...
    'Label', '  local | global', 'LabelOrientation', 'horizontal', 'FontSize', 9)

fprintf('\nDone. Figure generated.\n')

%% ---------------------------------------------------------------
%  LOCAL FUNCTIONS
% ----------------------------------------------------------------

function res = stiffness_eq47(n, ls_i, ls_v, Es, Is)
% STIFFNESS_EQ47  Evaluates Eq. (47) for finding end constraint factor n.
%
%  The equilibrium of rotations at the lattice nodes yields:
%    -2*s*c^2 + 2*s + s1/sqrt(2) + s1*c1/sqrt(2) = 0
%
%  s, c  : stability functions for inclined strut (COMPRESSION), Eq. (41)
%  s1, c1: stability functions for vertical strut (TENSION),     Eq. (43)
%
%  At buckling of the inclined strut: delta = pi*n (Eq. 42)
%  The vertical strut carries tension = sqrt(2) times the compression load.

    % --- Inclined strut (compression): stability functions (Eq. 41) ---
    delta  = pi*n;
    denom_s = 2 - 2*cos(delta) - delta*sin(delta);
    s       = delta*(sin(delta) - delta*cos(delta)) / denom_s;
    denom_c = sin(delta) - delta*cos(delta);
    c       = (delta - sin(delta)) / denom_c;

    % --- Vertical strut (tension): stability functions (Eq. 43) ---
    % Tension load on vertical strut = sqrt(2) * compression on inclined strut
    % At inclined strut buckling: F_inclined = n^2 * pi^2 * EI / ls_i^2
    FE_i   = pi^2 * Es * Is / ls_i^2;    % Euler load of inclined strut (n=1)
    F_v    = sqrt(2) * n^2 * FE_i;       % Tension in vertical strut [N]
    delta1 = sqrt(F_v/(Es*Is)) * ls_v;   % Non-dim. parameter (Eq. 44)

    denom_s1 = 2 - 2*cosh(delta1) + delta1*sinh(delta1);
    s1       = delta1*(delta1*cosh(delta1) - sinh(delta1)) / denom_s1;
    denom_c1 = delta1*cosh(delta1) - sinh(delta1);
    c1       = (sinh(delta1) - delta1) / denom_c1;

    % --- Eq. (47): nontrivial equilibrium condition ---
    res = -2*s*c^2 + 2*s + s1/sqrt(2) + s1*c1/sqrt(2);
end

function rgb = hex2rgb(hex)
% HEX2RGB  Converts '#RRGGBB' hex color string to [r g b] (range 0-1).
    hex = strrep(hex, '#', '');
    rgb = double(reshape(sscanf(hex, '%2x'), 3, []).') / 255;
end