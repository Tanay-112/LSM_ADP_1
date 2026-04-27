function core = unit_cell_FBCC(E0, d_strut, a_cell)
% FBCC core properties
%
% Inputs:
%   E0: young's modulus of the base material AlSi10Mg [N/mm^2]
%   d_strut: strut diameter [mm]
%   a_cell: cell size [mm]
%
% Outputs:
%   properties

    gamma = d_strut / a_cell;

    core.gamma = gamma;

    core.Ex = E0 * (1.661 * gamma^2.045);
    core.Ey = core.Ex;
    core.Ez = E0 * (2.442 * gamma^2.953);

    core.Gxy = E0 * (0.6748 * gamma^2.045);
    core.Gyz = E0 * (1.275  * gamma^2.015);
    core.Gzx = core.Gyz;

    core.vxy = 2.46 * gamma^2 + 1.526 * gamma - 0.5568;
    core.vyx = core.vxy;

    core.vxz = -2.551 * gamma^2 - 1.223 * gamma + 1.038;
    core.vyz = core.vxz;

    core.vzx = -1.27 * gamma^2 - 0.01504 * gamma + 0.6685;
    core.vzy = core.vzx;
end