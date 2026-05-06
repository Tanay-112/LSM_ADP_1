%% ========================================================================
%  3D FBCC LATTICE SANDWICH COLUMN - STRUT BUCKLING WITH SS-SS ENDS
%  ------------------------------------------------------------------------
%  This code is the 3D-FBCC conversion of a 2D-BCC strut-buckling idea.
%
%  Boundary condition used here:
%       Global sandwich column = simply supported at BOTH ends (SS-SS)
%
%  Important:
%  - This is a MATLAB analytical/frame proof model.
%  - Facesheets are represented by simplified beam-grillage restraints.
%  - For final publication-level validation, compare with Abaqus/ANSYS:
%       solid facesheets + beam lattice core + tie/coupling + eigen-buckling.
%
%  Main strut-buckling equation implemented for every selected candidate:
%
%       [Ke_i^G + Krest_i^3D] phi_i = Ncr_i [Kg_i^G] phi_i
%
%       Krest_i^3D = Krr - Krc * inv(Kcc) * Kcr
%
%       Pcr_i = Ncr_i / alpha_i
%
%  DOF order at each beam node:
%       [ux uy uz rx ry rz]
%
%  Element DOF order for one strut:
%       [uA vA wA rxA ryA rzA uB vB wB rxB ryB rzB]
%
%  Units:
%       mm, N, MPa = N/mm^2
% ========================================================================
function FBCC_3D_SS_SS_strut_buckling_full()
clc; close all;

%% ========================================================================
%  1. INPUTS: GEOMETRY, MATERIAL, AND NUMERICAL SETTINGS
% ========================================================================

% ---- AlSi10Mg / aluminium input ----
E0      = 70000;       % MPa = N/mm^2, solid AlSi10Mg Young's modulus input
nu0     = 0.33;        % Poisson's ratio input
G0      = E0/(2*(1+nu0));
sigmaY  = 250;         % MPa, optional yield cap only; edit if your data differs

% ---- FBCC unit cell and sandwich-column dimensions ----
a       = 4.0;         % mm, FBCC unit cell size
strut_d = 0.4;         % mm, circular strut diameter
Lcol    = 240.0;       % mm, sandwich column length along x
Bcol    = 24.0;        % mm, sandwich column width along y
Tcore   = 20.0;        % mm, lattice core thickness along z

% Number of cells. These must be integers for the generated lattice.
Nx = round(Lcol/a);
Ny = round(Bcol/a);
Nz = round(Tcore/a);

if abs(Nx*a - Lcol) > 1e-9 || abs(Ny*a - Bcol) > 1e-9 || abs(Nz*a - Tcore) > 1e-9
    error('Lcol, Bcol, and Tcore must be exact multiples of unit cell size a.');
end

% ---- simplified facesheet beam-grillage restraint ----
includeFaceGrillage = true;
tf = 2.0;             % mm, editable facesheet thickness for grillage restraint
faceStripWidth = a/2; % mm, spacing width assigned to each grillage beam

% ---- unit load for finding alpha_i ----
Punit = 1.0;           % N, total axial compression applied at right end

% ---- candidate control ----
numCandidates = 20;    % number of likely critical compressed struts checked by Krest
alphaTolerance = 1e-10;

% ---- support implementation ----
% Penalty springs allow every strut end to retain all 12 DOFs for Krest condensation.
% This makes the global SS-SS supports part of Krest.
supportPenalty = 1e13; % N/mm translational support spring; edit if numerical warnings appear

% ---- output filenames ----
outCSV = 'FBCC_3D_SS_SS_Krest_results.csv';
outMAT = 'FBCC_3D_SS_SS_Krest_workspace.mat';

%% ========================================================================
%  2. SECTION PROPERTIES AND FBCC EQUIVALENT PROPERTIES
% ========================================================================
A_core = pi*strut_d^2/4;
I_core = pi*strut_d^4/64;
J_core = pi*strut_d^4/32;
Nyield = A_core*sigmaY;

gamma = strut_d/a;
Ex_FBCC  = E0*(1.661*gamma^2.045);
Ey_FBCC  = Ex_FBCC;
Ez_FBCC  = E0*(2.442*gamma^2.953);
Gxy_FBCC = E0*(0.6748*gamma^2.045);
Gyz_FBCC = E0*(1.275*gamma^2.015);
Gzx_FBCC = Gyz_FBCC;
nuxy_FBCC = 2.46*gamma^2 + 1.526*gamma - 0.5568;
nuxz_FBCC = -2.551*gamma^2 - 1.223*gamma + 1.038;
nuzx_FBCC = -1.27*gamma^2 - 0.01504*gamma + 0.6685;

% Facesheet grillage section approximation.
% This is not a C3D20R solid facesheet. It is a simplified analytical restraint.
A_face = faceStripWidth*tf;
I_face = faceStripWidth*tf^3/12;
J_face = 2*I_face;

%% ========================================================================
%  3. GENERATE 3D FBCC LATTICE GEOMETRY
% ========================================================================
% FBCC unit cell used here:
%   - one body center per cell
%   - six face centers per cell
%   - eight corners per cell
%   - body center connected to 6 face centers  -> length a/2 = 2 mm
%   - body center connected to 8 corners       -> length sqrt(3)*a/2 = 3.464 mm
%
% This is the 3D conversion from the 2D-BCC idea:
%   2D BCC: one cell in x-z plane, inclined struts + vertical strut.
%   3D FBCC: cell is expanded in y-direction and face-center/body-center
%            nodes are included, giving struts in full x-y-z space.

fprintf('Generating 3D FBCC lattice with SS-SS global supports...\n');

nodeMap = containers.Map('KeyType','char','ValueType','int32');
nodes = zeros(0,3);       % physical coordinates [x y z]
nodeHK = zeros(0,3);      % integer half-step keys [ix iy iz]

maxCoreElem = Nx*Ny*Nz*14;
maxFaceElem = 200000; % safe upper bound; arrays truncated later
maxElem = maxCoreElem + maxFaceElem;

eN1 = zeros(maxElem,1);
eN2 = zeros(maxElem,1);
eIsCore = false(maxElem,1);
eType = strings(maxElem,1);
eA  = zeros(maxElem,1);
eIy = zeros(maxElem,1);
eIz = zeros(maxElem,1);
eJ  = zeros(maxElem,1);
eE  = zeros(maxElem,1);
eG  = zeros(maxElem,1);
nelem = 0;

% Helper for adding an element
    function addElement(n1,n2,isCore,typeName,Ae,Iye,Ize,Je,Ee,Ge)
        if n1 == n2
            return;
        end
        nelem = nelem + 1;
        if nelem > length(eN1)
            error('Element array too small. Increase maxFaceElem.');
        end
        eN1(nelem) = n1;
        eN2(nelem) = n2;
        eIsCore(nelem) = isCore;
        eType(nelem) = string(typeName);
        eA(nelem) = Ae;
        eIy(nelem) = Iye;
        eIz(nelem) = Ize;
        eJ(nelem) = Je;
        eE(nelem) = Ee;
        eG(nelem) = Ge;
    end

% Cell generation
for ix = 0:Nx-1
    for iy = 0:Ny-1
        for iz = 0:Nz-1
            % integer half-step coordinates
            bodyHK = [2*ix+1, 2*iy+1, 2*iz+1];
            body = getNodeID(bodyHK);

            faceHK = [
                2*ix,   2*iy+1, 2*iz+1;
                2*ix+2, 2*iy+1, 2*iz+1;
                2*ix+1, 2*iy,   2*iz+1;
                2*ix+1, 2*iy+2, 2*iz+1;
                2*ix+1, 2*iy+1, 2*iz;
                2*ix+1, 2*iy+1, 2*iz+2];

            for p = 1:6
                nf = getNodeID(faceHK(p,:));
                addElement(body,nf,true,'body-face',A_core,I_core,I_core,J_core,E0,G0);
            end

            cornerHK = [
                2*ix,   2*iy,   2*iz;
                2*ix+2, 2*iy,   2*iz;
                2*ix,   2*iy+2, 2*iz;
                2*ix+2, 2*iy+2, 2*iz;
                2*ix,   2*iy,   2*iz+2;
                2*ix+2, 2*iy,   2*iz+2;
                2*ix,   2*iy+2, 2*iz+2;
                2*ix+2, 2*iy+2, 2*iz+2];

            for p = 1:8
                nc = getNodeID(cornerHK(p,:));
                addElement(body,nc,true,'body-corner',A_core,I_core,I_core,J_core,E0,G0);
            end
        end
    end
end

nCoreElements = nelem;

%% ========================================================================
%  4. ADD SIMPLIFIED TOP/BOTTOM FACESHEET GRILLAGE
% ========================================================================
if includeFaceGrillage
    fprintf('Adding simplified top/bottom facesheet grillage...\n');

    for ii = 1:size(nodes,1)
        hk = nodeHK(ii,:);
        isBottom = (hk(3) == 0);
        isTop    = (hk(3) == 2*Nz);

        if ~(isBottom || isTop)
            continue;
        end

        % connect +x half-step neighbor in same z-plane
        hkx = hk + [1 0 0];
        keyx = makeKey(hkx);
        if isKey(nodeMap,keyx)
            jj = nodeMap(keyx);
            addElement(ii,jj,false,'face-grillage-x',A_face,I_face,I_face,J_face,E0,G0);
        end

        % connect +y half-step neighbor in same z-plane
        hky = hk + [0 1 0];
        keyy = makeKey(hky);
        if isKey(nodeMap,keyy)
            jj = nodeMap(keyy);
            addElement(ii,jj,false,'face-grillage-y',A_face,I_face,I_face,J_face,E0,G0);
        end
    end
end

% Truncate element arrays
eN1 = eN1(1:nelem);
eN2 = eN2(1:nelem);
eIsCore = eIsCore(1:nelem);
eType = eType(1:nelem);
eA  = eA(1:nelem);
eIy = eIy(1:nelem);
eIz = eIz(1:nelem);
eJ  = eJ(1:nelem);
eE  = eE(1:nelem);
eG  = eG(1:nelem);

nnode = size(nodes,1);
ndof = 6*nnode;

%% ========================================================================
%  5. ASSEMBLE GLOBAL 3D FRAME STIFFNESS MATRIX
% ========================================================================
fprintf('Assembling global 3D frame stiffness: %d nodes, %d elements, %d DOF...\n', nnode, nelem, ndof);

% each 12x12 element gives 144 entries
Iind = zeros(nelem*144,1);
Jind = zeros(nelem*144,1);
Vind = zeros(nelem*144,1);
ptr = 0;

% Store core element transformation and stiffness only when required later.
for e = 1:nelem
    n1 = eN1(e); n2 = eN2(e);
    x1 = nodes(n1,:); x2 = nodes(n2,:);
    Le = norm(x2-x1);
    if Le < 1e-12
        error('Zero-length element detected at e=%d.', e);
    end

    Keloc = beam3D_elastic_local(eE(e),eG(e),eA(e),eIy(e),eIz(e),eJ(e),Le);
    T = beam3D_transform(x1,x2);
    KeG = T.'*Keloc*T;
    edofs = elementDOF(n1,n2);

    [rr,cc] = ndgrid(edofs,edofs);
    rng = ptr + (1:144);
    Iind(rng) = rr(:);
    Jind(rng) = cc(:);
    Vind(rng) = KeG(:);
    ptr = ptr + 144;
end

Kglobal = sparse(Iind(1:ptr),Jind(1:ptr),Vind(1:ptr),ndof,ndof);
Kglobal = 0.5*(Kglobal + Kglobal.');

%% ========================================================================
%  6. GLOBAL SS-SS BOUNDARY CONDITIONS USING PENALTY SPRINGS
% ========================================================================
% Global sandwich simply supported at both ends for an axial column:
%   Left end x=0:
%       ux = uy = uz = 0  (reference/pin side; prevents rigid body motion)
%       rotations are free
%   Right end x=L:
%       uy = uz = 0
%       ux is free and receives axial compressive load
%       rotations are free
%
% This is a practical 3D-frame implementation of pinned-pinned / SS-SS ends.

leftNodes  = find(abs(nodes(:,1) - 0.0)  < 1e-9);
rightNodes = find(abs(nodes(:,1) - Lcol) < 1e-9);

supportDofs = [];
for id = leftNodes(:).'
    supportDofs = [supportDofs, dofOfNode(id,1), dofOfNode(id,2), dofOfNode(id,3)]; %#ok<AGROW>
end
for id = rightNodes(:).'
    supportDofs = [supportDofs, dofOfNode(id,2), dofOfNode(id,3)]; %#ok<AGROW>
end
supportDofs = unique(supportDofs);

Ksupport = sparse(supportDofs,supportDofs,supportPenalty*ones(length(supportDofs),1),ndof,ndof);
Kss = Kglobal + Ksupport;
Kss = 0.5*(Kss + Kss.');

% Unit compressive load at right end, distributed equally over all right nodes.
F = sparse(ndof,1);
for id = rightNodes(:).'
    F(dofOfNode(id,1)) = F(dofOfNode(id,1)) - Punit/length(rightNodes);
end

%% ========================================================================
%  7. STATIC SOLUTION TO FIND alpha_i FOR CORE STRUTS
% ========================================================================
fprintf('Solving unit-load static analysis for alpha_i...\n');
U = Kss \ F;

coreIDs = find(eIsCore);
nCore = length(coreIDs);
alpha = zeros(nCore,1);
Naxial = zeros(nCore,1);
Le_core = zeros(nCore,1);
NcrPin = zeros(nCore,1);
PscreenPin = inf(nCore,1);

for k = 1:nCore
    e = coreIDs(k);
    n1 = eN1(e); n2 = eN2(e);
    x1 = nodes(n1,:); x2 = nodes(n2,:);
    Le = norm(x2-x1);
    T = beam3D_transform(x1,x2);
    edofs = elementDOF(n1,n2);
    qloc = T*full(U(edofs));

    % Positive axial force = tension; negative = compression.
    Nax = eE(e)*eA(e)/Le * (qloc(7) - qloc(1));
    Ncomp = -Nax;

    Le_core(k) = Le;
    Naxial(k) = Nax;
    alpha(k) = Ncomp/Punit;

    NcrPin(k) = pi^2*E0*I_core/Le^2;
    if alpha(k) > alphaTolerance
        PscreenPin(k) = NcrPin(k)/alpha(k);
    end
end

compressedMask = alpha > alphaTolerance;
if ~any(compressedMask)
    error('No compressed struts found. Check loading direction or boundary conditions.');
end

[~,order] = sort(PscreenPin,'ascend');
order = order(isfinite(PscreenPin(order)));
numCandidates = min(numCandidates,length(order));
candidateCoreLocal = order(1:numCandidates);
candidateElemIDs = coreIDs(candidateCoreLocal);

%% ========================================================================
%  8. Krest CONDENSATION AND 12x12 STRUT EIGENVALUE BUCKLING
% ========================================================================
fprintf('Condensing Krest and solving 12x12 eigenproblem for %d candidate struts...\n', numCandidates);

results = table();

for c = 1:numCandidates
    e = candidateElemIDs(c);
    kLocal = candidateCoreLocal(c);
    n1 = eN1(e); n2 = eN2(e);
    x1 = nodes(n1,:); x2 = nodes(n2,:);
    Le = norm(x2-x1);
    edofs = elementDOF(n1,n2);

    % Target strut matrices in local and global coordinates.
    Keloc = beam3D_elastic_local(E0,G0,A_core,I_core,I_core,J_core,Le);
    Kgloc = beam3D_geometric_local_per_unit_force(Le);
    T = beam3D_transform(x1,x2);
    KeG = T.'*Keloc*T;
    KgG = T.'*Kgloc*T;

    % Remove target strut from remaining-sandwich stiffness.
    Krem = Kss;
    Krem(edofs,edofs) = Krem(edofs,edofs) - KeG;
    Krem = 0.5*(Krem + Krem.');

    % Static condensation onto the 12 target-strut DOFs.
    retained = edofs(:);
    allDofs = (1:ndof).';
    condensed = true(ndof,1);
    condensed(retained) = false;
    cdofs = allDofs(condensed);

    Krr = full(Krem(retained,retained));
    Krc = Krem(retained,cdofs);
    Kcr = Krem(cdofs,retained);
    Kcc = Krem(cdofs,cdofs);

    % Krest = Krr - Krc * inv(Kcc) * Kcr
    % Use Kcc\Kcr, not explicit inverse, for numerical stability.
    Krest = Krr - full(Krc * (Kcc \ Kcr));
    Krest = 0.5*(Krest + Krest.');

    % Generalized eigenvalue problem:
    %       (KeG + Krest) phi = Ncr KgG phi
    Aeig = 0.5*((KeG + Krest) + (KeG + Krest).');
    Beig = 0.5*(KgG + KgG.');

    lambda = eig(full(Aeig),full(Beig));
    lambda = real(lambda(abs(imag(lambda)) < 1e-6));
    lambda = lambda(isfinite(lambda));
    lambda = lambda(lambda > 1e-8);

    if isempty(lambda)
        NcrKrest = NaN;
        Keff = NaN;
        nFactor = NaN;
        Pcr = NaN;
    else
        NcrKrest = min(lambda);
        Keff = (pi/Le)*sqrt(E0*I_core/NcrKrest);
        nFactor = 1/Keff;
        Pcr = NcrKrest/alpha(kLocal);
    end

    Pyield = Nyield/alpha(kLocal);

    newRow = table(e,n1,n2,Le,alpha(kLocal),Naxial(kLocal),NcrPin(kLocal),PscreenPin(kLocal), ...
        NcrKrest,Keff,nFactor,Pcr,Nyield,Pyield, ...
        string(eType(e)), ...
        'VariableNames',{'strutID','nodeA','nodeB','L_mm','alpha_compression','Naxial_under_1N_N', ...
        'Ncr_pin_N','Pscreen_pin_N','Ncr_Krest_N','Keff','nFactor','Pcr_total_N', ...
        'Nyield_N','Pyield_total_N','strutType'});

    results = [results; newRow]; %#ok<AGROW>

    fprintf('Candidate %2d/%2d | eID=%6d | L=%8.5f mm | alpha=%10.4e | Ncr_Krest=%12.5f N | Keff=%8.4f | Pcr=%12.5f N\n', ...
        c,numCandidates,e,Le,alpha(kLocal),NcrKrest,Keff,Pcr);
end

% Pick final critical strut based on minimum total sandwich load Pcr.
valid = ~isnan(results.Pcr_total_N) & results.Pcr_total_N > 0;
if ~any(valid)
    error('No valid Krest eigenvalue result found. Try increasing supportPenalty or check geometry.');
end
[~,idxMin] = min(results.Pcr_total_N(valid));
validRows = find(valid);
criticalRow = results(validRows(idxMin),:);

%% ========================================================================
%  9. PRINT SUMMARY
% ========================================================================
fprintf('\n============================================================\n');
fprintf('INPUT GEOMETRY AND MATERIAL\n');
fprintf('============================================================\n');
fprintf('E0 AlSi10Mg input               = %.3f MPa\n',E0);
fprintf('nu0                             = %.5f\n',nu0);
fprintf('G0                              = %.3f MPa\n',G0);
fprintf('Unit cell a                     = %.3f mm\n',a);
fprintf('Strut diameter d                = %.3f mm\n',strut_d);
fprintf('gamma=d/a                       = %.6f\n',gamma);
fprintf('Column length L                 = %.3f mm (%d cells)\n',Lcol,Nx);
fprintf('Column width B                  = %.3f mm (%d cells)\n',Bcol,Ny);
fprintf('Core thickness Tc               = %.3f mm (%d cells)\n',Tcore,Nz);
fprintf('Nodes                           = %d\n',nnode);
fprintf('Core struts                     = %d\n',nCoreElements);
fprintf('All frame/grillage elements     = %d\n',nelem);
fprintf('Boundary condition              = SS-SS only\n');

fprintf('\n============================================================\n');
fprintf('CIRCULAR STRUT SECTION\n');
fprintf('============================================================\n');
fprintf('A                               = %.8f mm^2\n',A_core);
fprintf('I                               = %.8e mm^4\n',I_core);
fprintf('J                               = %.8e mm^4\n',J_core);
fprintf('A*sigmaY optional yield cap      = %.6f N\n',Nyield);

fprintf('\n============================================================\n');
fprintf('FBCC EQUIVALENT CORE PROPERTIES FROM PAPER TABLE FORMULAS\n');
fprintf('============================================================\n');
fprintf('Ex = Ey                         = %.6f MPa\n',Ex_FBCC);
fprintf('Ez                              = %.6f MPa\n',Ez_FBCC);
fprintf('Gxy                             = %.6f MPa\n',Gxy_FBCC);
fprintf('Gyz = Gzx                       = %.6f MPa\n',Gyz_FBCC);
fprintf('nuxy                            = %.6f\n',nuxy_FBCC);
fprintf('nuxz                            = %.6f\n',nuxz_FBCC);
fprintf('nuzx                            = %.6f\n',nuzx_FBCC);

fprintf('\n============================================================\n');
fprintf('STRUT LENGTH GROUPS - CORE ONLY\n');
fprintf('============================================================\n');
uniqueL = unique(round(Le_core,9));
for i = 1:length(uniqueL)
    Lval = uniqueL(i);
    count = sum(abs(Le_core - Lval) < 1e-7);
    NcrVal = pi^2*E0*I_core/Lval^2;
    fprintf('Length %.6f mm | count %6d | Euler pin-pin screening Ncr %.6f N\n',Lval,count,NcrVal);
end

fprintf('\n============================================================\n');
fprintf('FINAL CRITICAL STRUT - Krest BASED 12x12 EIGENVALUE RESULT\n');
fprintf('============================================================\n');
disp(criticalRow);

fprintf('Correct total-load formula used: Pcr_i = Ncr_i / alpha_i\n');
fprintf('Correct effective-length formula used: Keff = (pi/L)*sqrt(E*I/Ncr)\n');

%% ========================================================================
%  10. SAVE OUTPUTS AND PLOT CRITICAL STRUT
% ========================================================================
writetable(results,outCSV);
save(outMAT,'results','criticalRow','nodes','eN1','eN2','eIsCore','eType','E0','nu0','G0','a','strut_d','Lcol','Bcol','Tcore','Nx','Ny','Nz');

fprintf('\nSaved candidate result table to: %s\n',outCSV);
fprintf('Saved MATLAB workspace to:      %s\n',outMAT);

% Plot lattice skeleton and final critical strut.
figure('Color','w','Position',[80 80 900 550]);
hold on; axis equal; grid on; box on;
view(35,22);
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
title('3D FBCC lattice, SS-SS ends, critical strut highlighted');

% plot a reduced view for speed: all core struts in light grey
for e = coreIDs(:).'
    p1 = nodes(eN1(e),:); p2 = nodes(eN2(e),:);
    plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'-','Color',[0.75 0.75 0.75],'LineWidth',0.5);
end

ecrit = criticalRow.strutID;
p1 = nodes(eN1(ecrit),:); p2 = nodes(eN2(ecrit),:);
plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'r-','LineWidth',4);
scatter3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],80,'r','filled');
legend({'core struts','critical strut'},'Location','best');

fprintf('\nDONE.\n');

%% ========================================================================
%  LOCAL / NESTED FUNCTIONS
% ========================================================================

    function id = getNodeID(hk)
        key = makeKey(hk);
        if isKey(nodeMap,key)
            id = nodeMap(key);
        else
            id = int32(size(nodes,1) + 1);
            nodeMap(key) = id;
            nodes(double(id),:) = (a/2)*double(hk);
            nodeHK(double(id),:) = double(hk);
        end
    end

    function key = makeKey(hk)
        key = sprintf('%d_%d_%d',hk(1),hk(2),hk(3));
    end

    function dofs = elementDOF(n1,n2)
        dofs = [6*(n1-1)+(1:6), 6*(n2-1)+(1:6)];
    end

    function d = dofOfNode(n,localDof)
        d = 6*(n-1) + localDof;
    end

end % end main function

%% ========================================================================
%  STANDALONE LOCAL FUNCTIONS
% ========================================================================

function K = beam3D_elastic_local(E,G,A,Iy,Iz,J,L)
% Explicit 12x12 Euler-Bernoulli 3D beam stiffness matrix.
% Local DOF order:
% [u1 v1 w1 rx1 ry1 rz1 u2 v2 w2 rx2 ry2 rz2]

    a  = E*A/L;
    t  = G*J/L;

    by = 12*E*Iy/L^3;
    cy = 6*E*Iy/L^2;
    dy = 4*E*Iy/L;
    ey = 2*E*Iy/L;

    bz = 12*E*Iz/L^3;
    cz = 6*E*Iz/L^2;
    dz = 4*E*Iz/L;
    ez = 2*E*Iz/L;

    K = zeros(12,12);

    % axial u
    K(1,1)= a;  K(1,7)=-a;
    K(7,1)=-a;  K(7,7)= a;

    % torsion rx
    K(4,4)= t;  K(4,10)=-t;
    K(10,4)=-t; K(10,10)=t;

    % bending in local x-y plane: v and rz, using Iz
    ids = [2 6 8 12];
    K(ids,ids) = K(ids,ids) + [
         bz,  cz, -bz,  cz;
         cz,  dz, -cz,  ez;
        -bz, -cz,  bz, -cz;
         cz,  ez, -cz,  dz];

    % bending in local x-z plane: w and ry, using Iy
    ids = [3 5 9 11];
    K(ids,ids) = K(ids,ids) + [
         by, -cy, -by, -cy;
        -cy,  dy,  cy,  ey;
        -by,  cy,  by,  cy;
        -cy,  ey,  cy,  dy];

    K = 0.5*(K+K.');
end

function Kg = beam3D_geometric_local_per_unit_force(L)
% 12x12 geometric stiffness matrix per unit compressive axial force.
% Local DOF order:
% [u1 v1 w1 rx1 ry1 rz1 u2 v2 w2 rx2 ry2 rz2]
% Main contribution is flexural: v-rz plane and w-ry plane.

    Kg = zeros(12,12);

    kgPlane = (1/(30*L))*[
          36,    3*L,   -36,    3*L;
         3*L,  4*L^2,  -3*L,  -L^2;
         -36,   -3*L,    36,   -3*L;
         3*L,   -L^2,  -3*L,  4*L^2];

    % v-rz bending plane
    ids = [2 6 8 12];
    Kg(ids,ids) = Kg(ids,ids) + kgPlane;

    % w-ry bending plane: sign-adjusted to match elastic w-ry convention
    S = diag([1 -1 1 -1]);
    ids = [3 5 9 11];
    Kg(ids,ids) = Kg(ids,ids) + S*kgPlane*S;

    Kg = 0.5*(Kg+Kg.');
end

function T = beam3D_transform(x1,x2)
% Transformation q_local = T q_global.
% R rows are local axes expressed in global components.

    ex = (x2-x1)/norm(x2-x1);

    % choose reference vector not parallel to ex
    ref = [0 0 1];
    if abs(dot(ex,ref)) > 0.90
        ref = [0 1 0];
    end

    ey = cross(ref,ex);
    ey = ey/norm(ey);
    ez = cross(ex,ey);
    ez = ez/norm(ez);

    R = [ex; ey; ez];
    T = zeros(12,12);
    T(1:3,1:3)       = R;
    T(4:6,4:6)       = R;
    T(7:9,7:9)       = R;
    T(10:12,10:12)   = R;
end
