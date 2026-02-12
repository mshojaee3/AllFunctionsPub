function Tsummary = ShellResultantStress_Based_Volume_Integratin_of_RVE_10_2_26(csvFile, params, varargin)
% rveCsvToShellResultantStress_volume
% Build the same 3x8 summary table as your surface-based version, but using
% VOLUME averages (and volume-based slopes / moments) over the RVE.
%
% OUTPUT (table):
%   RowNames:  Input_prescribed, Input_calculated, Output
%   Cols:      11  22  12  31  23  11_1  22_1  12_1
%
% What it computes (volume-based):
%   - Fij_0 = volume-average of Fij over the RVE
%   - Fij_1 = best-fit slope wrt z using volume weighting:
%             K = ∫(f-f̄)(z-Z0)dV / ∫(z-Z0)^2 dV
%   - Pij_0 = volume-average of Pij over the RVE
%   - Output "_1" columns:
%       if useMomentOutputs=true (default):
%           MPij = (1/V) ∫ (z-Z0) Pij dV   (volume first moment)
%       else:
%           use volume-slope of P (Pij_1)
%
% Notes:
%   - Uses reference configuration volumes and reference z-centroids.
%   - Requires nodal CSV with X,Y,Z and U_U1,U_U2,U_U3; stresses S** for P.
%
% Options:
%   'tol'               (default 1e-8)
%   'Z0_center'        (default 0.5)
%   'useMomentOutputs'  true(default): Output "_1" are volume moments; false: use P slopes
%   'thrFactor'         (default 1e-6) row-wise tiny->0 threshold factor

% -------------------- parse inputs --------------------
p = inputParser;
p.addRequired('csvFile', @(s) ischar(s) || isstring(s));
p.addRequired('params', @(s) isstruct(s));

p.addParameter('tol', 1e-8, @(x) isscalar(x) && x>0);
p.addParameter('Z0_center', 0.5, @(x) isscalar(x) && isnumeric(x));
p.addParameter('useMomentOutputs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('thrFactor', 1e-6, @(x) isscalar(x) && x>0);
p.parse(csvFile, params, varargin{:});

tol       = p.Results.tol;
Z0_center= p.Results.Z0_center;
useMomentOutputs = p.Results.useMomentOutputs;
thrFactor = p.Results.thrFactor;

csvFile = char(csvFile);
assert(isfile(csvFile), 'CSV not found: %s', csvFile);

% Validate required params fields
req = {'H11_in','H22_in','H12_in','H31_in','H32_in','K11_in','K22_in','K12_in'};
for i=1:numel(req)
    assert(isfield(params,req{i}), 'params.%s is missing', req{i});
end

% -------------------- read csv --------------------
T = readtable(csvFile);

% -------------------- 1) Coordinates --------------------
X = T.X(:); Y = T.Y(:); Z = T.Z(:);
V = [X Y Z];

% -------------------- 2) Structured grid -> connectivity --------------------
ux = unique(round(X/tol)*tol); ux = sort(ux);
uy = unique(round(Y/tol)*tol); uy = sort(uy);
uz = unique(round(Z/tol)*tol); uz = sort(uz);

nx = numel(ux); ny = numel(uy); nz = numel(uz);
if nx*ny*nz ~= numel(X)
    warning(['nx*ny*nz ~= #nodes. Mesh might not be perfectly structured.\n' ...
        'Try increasing tol (e.g. 1e-6).']);
end

key3 = @(x,y,z) sprintf('%.12g_%.12g_%.12g', round(x/tol)*tol, round(y/tol)*tol, round(z/tol)*tol);
mp = containers.Map('KeyType','char','ValueType','double');
for r = 1:numel(X)
    mp(key3(X(r),Y(r),Z(r))) = r;
end

nodeId = zeros(nx, ny, nz);
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            kk = key3(ux(i), uy(j), uz(k));
            if ~isKey(mp, kk)
                error('Missing node at (%g,%g,%g). Increase tol or mesh is not structured.', ux(i),uy(j),uz(k));
            end
            nodeId(i,j,k) = mp(kk);
        end
    end
end

ne = (nx-1)*(ny-1)*(nz-1);
conn = zeros(ne, 8);
e = 0;
for i = 1:(nx-1)
    for j = 1:(ny-1)
        for k = 1:(nz-1)
            e = e + 1;
            n1 = nodeId(i,   j,   k);
            n2 = nodeId(i+1, j,   k);
            n3 = nodeId(i+1, j+1, k);
            n4 = nodeId(i,   j+1, k);
            n5 = nodeId(i,   j,   k+1);
            n6 = nodeId(i+1, j,   k+1);
            n7 = nodeId(i+1, j+1, k+1);
            n8 = nodeId(i,   j+1, k+1);
            conn(e,:) = [n1 n2 n3 n4 n5 n6 n7 n8];
        end
    end
end

% -------------------- element volumes (reference) + element centroids (reference) --------------------
Ve  = elementVolumesC3D8_ref(V, conn);
Vtot = sum(Ve);
Xe  = elementCentroidsC3D8_ref(V, conn);
ze  = Xe(:,3);

% -------------------- 5) Collect element-wise fields into "out" --------------------
out = struct();

fieldMap = { ...
    'S11',  'S11'; ...
    'S22',  'S22'; ...
    'S33',  'S33'; ...
    'S12',  'S12'; ...
    'S13',  'S13'; ...
    'S23',  'S23'; ...
    'LE11', 'LE11'; ...
    'LE22', 'LE22'; ...
    'LE33', 'LE33'; ...
    'LE12', 'LE12'; ...
    'LE13', 'LE13'; ...
    'LE23', 'LE23'; ...
    'U1',   'U_U1'; ...
    'U2',   'U_U2' ...
    };

for k = 1:size(fieldMap,1)
    label = fieldMap{k,1};
    col   = fieldMap{k,2};
    if ismember(col, T.Properties.VariableNames)
        Fn = T.(col)(:);
        out.(label) = mean(Fn(conn), 2); % nodal->element average
    end
end

% -------------------- 6) Compute F and detF (C3D8R centroid) --------------------
mustCols = {'U_U1','U_U2','U_U3'};
assert(all(ismember(mustCols, T.Properties.VariableNames)), ...
    'CSV is missing displacement columns U_U1/U_U2/U_U3.');

U1 = T.U_U1(:); U2 = T.U_U2(:); U3 = T.U_U3(:);
U  = [U1 U2 U3];

hex_dN = @(xi,eta,zeta) 0.125 * [ ...
    -(1-eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1-eta);  % N1
     (1-eta)*(1-zeta),  -(1+xi)*(1-zeta),  -(1+xi)*(1-eta);  % N2
     (1+eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1+eta);  % N3
    -(1+eta)*(1-zeta),   (1-xi)*(1-zeta),  -(1-xi)*(1+eta);  % N4
    -(1-eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1-eta);  % N5
     (1-eta)*(1+zeta),  -(1+xi)*(1+zeta),   (1+xi)*(1-eta);  % N6
     (1+eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1+eta);  % N7
    -(1+eta)*(1+zeta),   (1-xi)*(1+zeta),   (1-xi)*(1+eta)]; % N8
dN = hex_dN(0,0,0);

F    = zeros(ne, 9);
detF = zeros(ne, 1);
for ee = 1:ne
    ids = conn(ee,:);
    Xe_n  = V(ids,:);
    Ue  = U(ids,:);
    xe  = Xe_n + Ue;

    J0 = Xe_n.' * dN;
    J1 = xe.'   * dN;
    if rcond(J0) < 1e-14
        error('Element %d has near-singular J0. Check mesh/ordering.', ee);
    end
    Fe = J1 / J0;

    F(ee,:)  = [Fe(1,1) Fe(1,2) Fe(1,3) ...
                Fe(2,1) Fe(2,2) Fe(2,3) ...
                Fe(3,1) Fe(3,2) Fe(3,3)];
    detF(ee) = det(Fe);
end

out.F11 = F(:,1); out.F12 = F(:,2); out.F13 = F(:,3);
out.F21 = F(:,4); out.F22 = F(:,5); out.F23 = F(:,6);
out.F31 = F(:,7); out.F32 = F(:,8); out.F33 = F(:,9);
out.detF = detF;

% -------------------- 8) Compute P (1st PK) from Cauchy stress --------------------
needS = all(ismember({'S11','S22','S33','S12','S13','S23'}, T.Properties.VariableNames));
if needS
    s11n = T.S11(:); s22n = T.S22(:); s33n = T.S33(:);
    s12n = T.S12(:); s13n = T.S13(:); s23n = T.S23(:);

    s11e = mean(s11n(conn),2); s22e = mean(s22n(conn),2); s33e = mean(s33n(conn),2);
    s12e = mean(s12n(conn),2); s13e = mean(s13n(conn),2); s23e = mean(s23n(conn),2);

      % Preallocate ALL 9 components
      P11=zeros(ne,1); P12=zeros(ne,1); P13=zeros(ne,1);
      P21=zeros(ne,1); P22=zeros(ne,1); P23=zeros(ne,1);
      P31=zeros(ne,1); P32=zeros(ne,1); P33=zeros(ne,1);

    for ee = 1:ne
        Fe = [F(ee,1) F(ee,2) F(ee,3);
              F(ee,4) F(ee,5) F(ee,6);
              F(ee,7) F(ee,8) F(ee,9)];
        J = detF(ee);

        sigma = [s11e(ee) s12e(ee) s13e(ee);
                 s12e(ee) s22e(ee) s23e(ee);
                 s13e(ee) s23e(ee) s33e(ee)];

        Finv = inv(Fe);
        Pmat = J * sigma * (Finv.');   % P = J*sigma*F^{-T}

       % Store all components
        P11(ee)=Pmat(1,1); P12(ee)=Pmat(1,2); P13(ee)=Pmat(1,3);
        P21(ee)=Pmat(2,1); P22(ee)=Pmat(2,2); P23(ee)=Pmat(2,3);
        P31(ee)=Pmat(3,1); P32(ee)=Pmat(3,2); P33(ee)=Pmat(3,3);
    end

       % Output ALL 9
    out.P11 = P11; out.P12 = P12; out.P13 = P13;
    out.P21 = P21; out.P22 = P22; out.P23 = P23;
    out.P31 = P31; out.P32 = P32; out.P33 = P33;
else
    warning('Stress columns S** not found -> output P** will be NaN.');
    out.P11 = nan(ne,1); out.P12 = nan(ne,1); out.P13 = nan(ne,1);
    out.P21 = nan(ne,1); out.P22 = nan(ne,1); out.P23 = nan(ne,1);
    out.P31 = nan(ne,1); out.P32 = nan(ne,1); out.P33 = nan(ne,1);
end

% ===================== VOLUME extraction (F0,F1,P0,P1) =====================
% ---- 0th order: volume averages ----
% ---- 0th order: volume averages (ALL components) ----
% F_0 (9 components)
F11_0 = volumeAvg('F11');  F12_0 = volumeAvg('F12');  F13_0 = volumeAvg('F13');
F21_0 = volumeAvg('F21');  F22_0 = volumeAvg('F22');  F23_0 = volumeAvg('F23');
F31_0 = volumeAvg('F31');  F32_0 = volumeAvg('F32');  F33_0 = volumeAvg('F33');

% P_0 (9 components)
P11_0 = volumeAvg('P11');  P12_0 = volumeAvg('P12');  P13_0 = volumeAvg('P13');
P21_0 = volumeAvg('P21');  P22_0 = volumeAvg('P22');  P23_0 = volumeAvg('P23');
P31_0 = volumeAvg('P31');  P32_0 = volumeAvg('P32');  P33_0 = volumeAvg('P33');

% ---- 1st order: volume-based slopes wrt z (best-fit) ----
F11_1 = volumeSlope('F11'); % K11-like
F22_1 = volumeSlope('F22'); % K22-like
F12_1 = volumeSlope('F12'); % K12-like

P11_1 = volumeSlope('P11');
P22_1 = volumeSlope('P22');
P12_1 = volumeSlope('P12');

% ---- Optional: volume "moment" of P (first moment about Z0_center) ----
MP11 = NaN; MP22 = NaN; MP12 = NaN;
if useMomentOutputs && needS
    MP11 = volumeFirstMoment('P11');
    MP22 = volumeFirstMoment('P22');
    MP12 = volumeFirstMoment('P12');
end

% ===================== Build 3x8 summary table =====================
row_input = [params.H11_in, params.H22_in, params.H12_in, params.H31_in, params.H32_in, ...
             params.K11_in, params.K22_in, params.K12_in];

row_calc_input = [F11_0, F22_0, F12_0, F31_0, F32_0, F11_1, F22_1, F12_1];

if useMomentOutputs
    row_output = [P11_0, P22_0, P12_0, P31_0, P32_0, MP11, MP22, MP12];
else
    row_output = [P11_0, P22_0, P12_0, P31_0, P32_0, P11_1, P22_1, P12_1];
end

colNames = {'11','22','12','31','32','11_1','22_1','12_1'};
rowNames = {'Input_prescribed','Input_calculated','Output'};
M = [row_input; row_calc_input; row_output];

% ---- row-wise tiny->0 ----
for r = 1:size(M,1)
    big = max(abs(M(r,:)));
    if big > 0
        thr = thrFactor * big;
        M(r, abs(M(r,:)) < thr) = 0;
    end
end

Tsummary = array2table(M, 'VariableNames', colNames, 'RowNames', rowNames);

% disp('===== 3x8 SUMMARY (VOLUME-BASED, ROW-WISE SMALL->0) =====');
% disp(Tsummary);

% -------------------- nested volume helpers --------------------
    function avg = volumeAvg(fieldName)
        fE = out.(fieldName)(:);
        avg = sum(fE .* Ve) / max(Vtot, eps);
    end

    function K = volumeSlope(fieldName)
        fE = out.(fieldName)(:);
        Havg = sum(fE .* Ve) / max(Vtot, eps);
        zrel = (ze - Z0_center);
        num = sum( (fE - Havg) .* zrel .* Ve );
        den = sum( (zrel.^2) .* Ve );
        K = num / max(den, eps);
    end

    function Mv = volumeFirstMoment(fieldName)
        fE = out.(fieldName)(:);
        zrel = (ze - Z0_center);
        Mv = sum( zrel .* fE .* Ve ) / max(Vtot, eps);
    end

end

% =====================================================================
% Reference element volumes: C3D8 with 2x2x2 Gauss quadrature
% =====================================================================
function Ve = elementVolumesC3D8_ref(V, conn)
ne = size(conn,1);
Ve = zeros(ne,1);

g = 1/sqrt(3);
gp = [-g, +g];

hex_dN = @(xi,eta,zeta) 0.125 * [ ...
    -(1-eta)*(1-zeta),  -(1-xi)*(1-zeta),  -(1-xi)*(1-eta);  % N1
     (1-eta)*(1-zeta),  -(1+xi)*(1-zeta),  -(1+xi)*(1-eta);  % N2
     (1+eta)*(1-zeta),   (1+xi)*(1-zeta),  -(1+xi)*(1+eta);  % N3
    -(1+eta)*(1-zeta),   (1-xi)*(1-zeta),  -(1-xi)*(1+eta);  % N4
    -(1-eta)*(1+zeta),  -(1-xi)*(1+zeta),   (1-xi)*(1-eta);  % N5
     (1-eta)*(1+zeta),  -(1+xi)*(1+zeta),   (1+xi)*(1-eta);  % N6
     (1+eta)*(1+zeta),   (1+xi)*(1+zeta),   (1+xi)*(1+eta);  % N7
    -(1+eta)*(1+zeta),   (1-xi)*(1+zeta),   (1-xi)*(1+eta)]; % N8

for ee = 1:ne
    ids = conn(ee,:);
    Xe  = V(ids,:);  % 8x3 reference nodes
    vol = 0;
    for i = 1:2
        xi = gp(i);
        for j = 1:2
            eta = gp(j);
            for k = 1:2
                zeta = gp(k);
                dN = hex_dN(xi,eta,zeta);   % 8x3
                J  = Xe.' * dN;             % 3x3
                vol = vol + det(J);         % weights are 1
            end
        end
    end
    Ve(ee) = vol;
end
end

% =====================================================================
% Reference element centroids (simple nodal average)
% =====================================================================
function Xe = elementCentroidsC3D8_ref(V, conn)
ne = size(conn,1);
Xe = zeros(ne,3);
for ee = 1:ne
    Xe(ee,:) = mean(V(conn(ee,:),:), 1);
end
end



