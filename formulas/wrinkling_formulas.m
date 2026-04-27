function [Pcr_wrinkling, details] = wrinkling_formulas(params)
% wrinkling_formulas use Plantema's Formula
%
% Inputs:
%   C : constant, Plantema C = 1.7 is used by default
%   b: width of column [mm]
%   f: thickness of facesheet [mm]
%   Ef: young's modulus of facesheet [N/mm^2]
%
% Outputs:
%   Pcr_wrinkling  - wrinkling buckling load [N]
%   details - struct containing intermediate quantities

%C  = params.C; 
C  = 1.7;  % Plantema coefficient
b  = params.b;
f  = params.f;
Ef = params.Ef;

% get the FBCC core properties
core = unit_cell_FBCC(params.E0, params.d_strut, params.a_cell);
Ec_z = core.Ez;
Gc   = core.Gyz;

Pcr_w_per_width = C * f *(Ef * Ec_z * Gc)^(1/3);
Pcr_wrinkling = Pcr_w_per_width * b;

details.Pcr_w_per_width = Pcr_w_per_width;

details.Ec_z = Ec_z;
details.Gc   = Gc;

details.C = C;
details.core = core;

end