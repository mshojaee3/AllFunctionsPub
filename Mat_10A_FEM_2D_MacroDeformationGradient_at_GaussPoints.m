function [coor, Fij, Labels_np, out] = Mat_10A_FEM_2D_MacroDeformationGradient_at_GaussPoints(parms, RUC_Params)

%RUN_FULLFIELD_MACROF  Full-field Abaqus -> compute macro deformation gradient F at Gauss points
%
% Inputs
%   parms       : struct (must include .load_case, .JOB or , .NX,.NY,.Lx,.Ly,.mesh_size, etc.)
%   RUC_Params  : struct (must include .run.simDir, .run.pyFile, .Lx,.Ly,.mesh_size, .BC_MODE, etc.)
%
% Outputs
%   coor      : [nGP x 2] global Gauss point coords = [x_gp, y_gp]
%   Fij       : cell {F11,F12,F21,F22}, each [nGP x 1]
%   Labels_np : [nNode x 1] nodal labels read from CSV
%   out       : struct with extra internals (Hmacro, Grad_U1_gp, Grad_U2_gp, L, Results, idxLgp)

res_name = sprintf('RUC_results.mat');
run_ruc = true;
simDir = parms.simDir;
master_py=parms.master_py;
abaqus_cmd=parms.abaqus_cmd;
Lx_ruc    = RUC_Params.Lx;
Ly_ruc    = RUC_Params.Ly;
mesh_size = RUC_Params.mesh_size;

if run_ruc


    [L, Results] = AllFunctions.Mat_9A_FEM_2D_LocalizationTensorFromCorrectors_at_GaussPoints(RUC_Params);

    Jobs = {'Uniaxialx', 'Uniaxialy', 'Shear12', 'Shear21'};

    for k = 1:4
        job_name = Jobs{k};

 
        grads = {Results(k).G11, Results(k).G12, Results(k).G21, Results(k).G22};

    end

    save(res_name, 'L', 'Results');
    save('LocalizationTensor_RUC.mat','L');
end

%% -------------------------------------------------------------------------
%% 3. FULL-FIELD SIMULATION (Macroscopic)
%% -------------------------------------------------------------------------

AllFunctions.Mat_5A_safeCleanRunDir(simDir);

% 2.1 Generate Python Script using YOUR function
out_py_name = [parms.JOB, '_run.py'];
full_py_path = fullfile(simDir, out_py_name);

fprintf('>> Generating script using AllFunctions.Mat_7A_updatePyFromParams...\n');
AllFunctions.Mat_7A_updatePyFromParams(parms, master_py, full_py_path);


% 2.2 Run Abaqus
cmd = sprintf('cmd /c "cd /d ""%s"" && %s cae noGUI=""%s"""', ...
    fullfile(pwd, simDir), abaqus_cmd, out_py_name);
status = system(cmd);

if status ~= 0
    error('Abaqus simulation failed. Check .log files in %s', simDir);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
csv_node  = fullfile(simDir, sprintf('%s_NodalData.csv', parms.JOB));
csv_conn  = fullfile(simDir, sprintf('%s_Connectivity.csv', parms.JOB));
nodal_table = readtable(csv_node);

raw_conn = readmatrix(csv_conn, 'NumHeaderLines', 1);
Connectivity = raw_conn;

Labels_np = nodal_table.Node_Label;
X_np      = nodal_table.X;
Y_np      = nodal_table.Y;
U1_np     = nodal_table.U_U1;
U2_np   = nodal_table.U_U2;

% 1. Calculate Gradients for U1 (Displacement X)
Grad_U1_gp = AllFunctions.Mat_2C_FEM_2D_GaussPoint_gradient_from_NodalPoint(X_np, Y_np, U1_np, Connectivity, Labels_np);

% 2. Calculate Gradients for U2 (Displacement Y)
Grad_U2_gp = AllFunctions.Mat_2C_FEM_2D_GaussPoint_gradient_from_NodalPoint(X_np, Y_np, U2_np, Connectivity, Labels_np);

% 3. Extract Gradient Components for F Matrix
% H11 = dU1/dX, H12 = dU1/dY
H11_gp_ab = Grad_U1_gp(:, 6); 
H12_gp_ab = Grad_U1_gp(:, 7);

% H21 = dU2/dX, H22 = dU2/dY
H21_gp_ab = Grad_U2_gp(:, 6);
H22_gp_ab = Grad_U2_gp(:, 7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% -------------------------------------------------------------------------
%% 4. Macrostrain calculation (homogenized strain) 
%% -------------------------------------------------------------------------


%% --- Map global GP -> local RUC coordinates ---

y1loc = mod(Grad_U1_gp(:, 3), Lx_ruc);
y2loc = mod(Grad_U1_gp(:, 4), Ly_ruc);

tolKey = max(1e-8, 1e-4*mesh_size);

% Snap upper boundary back to zero
y1loc(abs(y1loc - Lx_ruc) < tolKey) = 0;
y2loc(abs(y2loc - Ly_ruc) < tolKey) = 0;

keyInt = @(a,b) strcat(string(round(a./tolKey)), "_", string(round(b./tolKey)));

% Build RUC key map
keys_ruc = keyInt(L.Y1, L.Y2);
[uniq_keys, ia] = unique(keys_ruc, 'stable');
idxMap = containers.Map(uniq_keys, ia);

% Find corresponding RUC index for each GP
keys_gp = keyInt(y1loc, y2loc);
x_gp = Grad_U1_gp(:,3);
y_gp = Grad_U1_gp(:,4);
nGP  = size(Grad_U1_gp,1);
idxLgp = zeros(nGP,1);
missing = false(nGP,1);
coor = [x_gp, y_gp];
for i = 1:nGP
    k = keys_gp(i);
    if isKey(idxMap, k)
        idxLgp(i) = idxMap(k);
    else
        missing(i) = true;
        [~, idxLgp(i)] = min((L.Y1 - y1loc(i)).^2 + (L.Y2 - y2loc(i)).^2);
    end
end

fprintf('   -> GP mapping: %d/%d used nearest fallback.\n', nnz(missing), nGP);

%% --- Pull localization tensor at mapped points ---
L11E11 = L.L11_E11(idxLgp);  L11E22 = L.L11_E22(idxLgp);  L11E12 = L.L11_E12(idxLgp); L11E21 = L.L11_E21(idxLgp);
L12E11 = L.L12_E11(idxLgp);  L12E22 = L.L12_E22(idxLgp);  L12E12 = L.L12_E12(idxLgp); L12E21 = L.L12_E21(idxLgp);
L21E11 = L.L21_E11(idxLgp);  L21E22 = L.L21_E22(idxLgp);  L21E12 = L.L21_E12(idxLgp); L21E21 = L.L21_E21(idxLgp);
L22E11 = L.L22_E11(idxLgp);  L22E22 = L.L22_E22(idxLgp);  L22E12 = L.L22_E12(idxLgp); L22E21 = L.L22_E21(idxLgp);


%% -------------------------------------------------------------------------
%% 4.2  Local regression to compute H_macro from  H_full = L * H_macro
%%      (macro varies linearly inside each neighborhood)
%% -------------------------------------------------------------------------

% ---- Local weighted LS settings ----
kNN    = 60;     % neighborhood size (>= ~40 recommended since we solve 12 unknowns)
lambda = 1e-8;   % ridge regularization (try 1e-10 ... 1e-6)
p      = 2;      % weight exponent (used with distance, not squared distance)

Hmacro11 = zeros(nGP,1);
Hmacro12 = zeros(nGP,1);
Hmacro21 = zeros(nGP,1);
Hmacro22 = zeros(nGP,1);

% Precompute squared distances between GPs
XY = [x_gp(:), y_gp(:)];
D2 = pdist2(XY, XY, 'squaredeuclidean');

for g = 1:nGP

    % --- k nearest neighbors of point g (including itself) ---
    [d2_sorted, idx_sorted] = sort(D2(g,:), 'ascend');
    nb = idx_sorted(1:min(kNN, nGP));
    d2 = d2_sorted(1:numel(nb));      % squared distances

    % --- weights: use distance-based kernel (Gaussian-like) ---
    d = sqrt(d2);                      % convert to distances
    h = median(d(2:end));              % scale (ignore zero distance to itself)
    if ~isfinite(h) || h <= 0
        h = max(d(2), eps);
    end
    w = exp(-(d./h).^p);
    Wsqrt = sqrt(w(:));

    % --- build local system for theta (12 unknowns) ---
    % theta = [a11 b11 c11 a12 b12 c12 a21 b21 c21 a22 b22 c22]^T
    % Hij_macro(x,y) ~ aij + bij*(x-xg) + cij*(y-yg)
    m = numel(nb);
    A_loc = zeros(4*m, 12);
    b_loc = zeros(4*m, 1);

    x0 = x_gp(g);
    y0 = y_gp(g);

    for t = 1:m
        hIdx = nb(t);

        % A(h) = L(h) (4x4), b(h)=H_full(h) (4x1)
        Ah = [ L11E11(hIdx), L11E12(hIdx), L11E21(hIdx), L11E22(hIdx);
               L12E11(hIdx), L12E12(hIdx), L12E21(hIdx), L12E22(hIdx);
               L21E11(hIdx), L21E12(hIdx), L21E21(hIdx), L21E22(hIdx);
               L22E11(hIdx), L22E12(hIdx), L22E21(hIdx), L22E22(hIdx) ];

        bh = [ H11_gp_ab(hIdx);
               H12_gp_ab(hIdx);
               H21_gp_ab(hIdx);
               H22_gp_ab(hIdx) ];

        dx = x_gp(hIdx) - x0;
        dy = y_gp(hIdx) - y0;
        Ph = [1, dx, dy];   % 1x3

        % Mh maps theta -> x(h) = [H11 H12 H21 H22]^T at neighbor h (4x12)
        Mh = [ Ph, zeros(1,9);
               zeros(1,3), Ph, zeros(1,6);
               zeros(1,6), Ph, zeros(1,3);
               zeros(1,9), Ph ];

        % equation at neighbor h:  b(h) = A(h)*x(h) = A(h)*Mh*theta
        row = (t-1)*4 + (1:4);
        A_loc(row,:) = Ah * Mh;
        b_loc(row)   = bh;
    end

    % --- apply weights (scale each 4-row block by sqrt(w_t)) ---
    for t = 1:m
        row = (t-1)*4 + (1:4);
        A_loc(row,:) = Wsqrt(t) * A_loc(row,:);
        b_loc(row)   = Wsqrt(t) * b_loc(row);
    end

    % --- ridge solve for theta ---
    theta = (A_loc.'*A_loc + lambda*eye(12)) \ (A_loc.'*b_loc);

    % --- H_macro at center point g (dx=dy=0) is the constant terms aij ---
    Hmacro11(g) = theta(1);
    Hmacro12(g) = theta(4);
    Hmacro21(g) = theta(7);
    Hmacro22(g) = theta(10);
end

% -----------------------------
% 7) Build F = I + Hmacro
% -----------------------------
F11 = 1 + Hmacro11;
F12 = Hmacro12;
F21 = Hmacro21;
F22 = 1 + Hmacro22;

Fij = {F11, F12, F21, F22};


% Optional extra outputs
out = struct();
out.test = 1;

