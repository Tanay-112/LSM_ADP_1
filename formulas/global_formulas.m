function [Pcr_global, mode_global, details] = global_formulas(params)
% global_formulas use Allen's Formula - thin mode
%
% Inputs:
%   b: width of column [mm]
%   c: thickness of core [mm]
%   d: face-sheet mid-surface spacing [mm]
%   f: thickness of facesheet [mm]
%   Ef: young's modulus of facesheet [N/mm^2]
%   K: coefficient of BCs
%   L: length of column [mm]
%
% Outputs:
%   Pcr_global  - governing global buckling load [N]
%   mode_global - governing global buckling mode
%   details - struct containing intermediate quantities

b  = params.b;
c  = params.c;
d  = params.d;
f  = params.f;
Ef = params.Ef;
K  = params.K;
L  = params.L;

% get the FBCC core properties
core = unit_cell_FBCC(params.E0, params.d_strut, params.a_cell);
Gc_out = core.Gzx;
Gc_in  = core.Gxy;
Ec_i   = core.Ex;

A = b * d^2 / c;

%   Ds_out: out-of-plane buckling
Ds_out = Ef * b * f *d^2 / 2;

%   Ds_in: in-plane buckling
Ds_in_paper = Ef * f * b^3 / 6 + Ec_i * b^3 / 12;
Ds_in_dim   = Ef * f * b^3 / 6 + Ec_i * c * b^3 / 12;
Ds_in = Ds_in_dim;
% choose dimensionally consistent version

%   PE: Euler-type global buckling load [N]
PE_out = pi^2 * Ds_out / (K * L)^2;
PE_in = pi^2 * Ds_in / (K * L)^2;

Pcr_out = A * Gc_out * PE_out / (A * Gc_out + PE_out);
Pcr_in  = A * Gc_in  * PE_in  / (A * Gc_in  + PE_in);

if Pcr_out <= Pcr_in
        Pcr_global = Pcr_out;
        mode_global = "global_out_of_plane";
    else
        Pcr_global = Pcr_in;
        mode_global = "global_in_plane";
end

details.Pcr_out = Pcr_out;
details.Pcr_in  = Pcr_in;
details.PE_out  = PE_out;
details.PE_in   = PE_in;

details.Ds_out = Ds_out;

details.Ds_in_paper = Ds_in_paper;
details.Ds_in_dim   = Ds_in_dim;
details.Ds_in       = Ds_in;

details.Gc_out = Gc_out;
details.Gc_in  = Gc_in;
details.Ec_i   = Ec_i;

details.core  = core;
details.gamma = core.gamma;

end