function [L, Results] = Mat_9A_FEM_2D_LocalizationTensorFromCorrectors_at_GaussPoints(params, pyFile, abaqus_cmd, simDir, res_name, doPlot)
% Mat_9A_FEM_2D_LocalizationTensorFromCorrectors_at_GaussPoints
% Runs 4 RUC corrector problems (11,22,12,21), computes corrector gradients
% at Gauss points, and assembles the localization tensor L at Gauss points.
%
% INPUTS
%   params      : struct with (at least) fields:
%                 params.Lx, params.Ly, params.mesh_size)
%   pyFile      : master Abaqus python file (e.g., 'Main2D.py')
%   abaqus_cmd  : Abaqus command (e.g., 'abaqus')
%   simDir      : working folder (e.g., 'WorkDir_Abaqus')
%   res_name    : .mat filename to save L (e.g., 'RUC_results.mat') (optional)
%   doPlot      : true/false for plots (optional, default false)
%
% OUTPUTS
%   L       : struct with fields Y1,Y2,NodeLabels and Lij_Ekl components
%   Results : struct array (1..4) with Gauss point coords and Gij correctors
%
% DEPENDS ON (your package functions)
%   AllFunctions.Mat_5A_safeCleanRunDir
%   AllFunctions.Mat_7A_updatePyFromParams
%   AllFunctions.Mat_2C_FEM_2D_GaussPoint_gradient_from_NodalPoint

% ----------------- defaults / validate -----------------
assert(isstruct(params), 'params must be a struct.');
assert(isfield(params,'Lx') && isfield(params,'Ly') && isfield(params,'mesh_size'), ...
    'params must contain: Lx, Ly, mesh_size');
    
PBC_VALUE = 'ON';

assert(ischar(pyFile) || isstring(pyFile), 'pyFile must be char/string.');
pyFile = char(pyFile);
assert(isfile(pyFile), 'Cannot find python file: %s', pyFile);

if nargin < 3 || isempty(abaqus_cmd), abaqus_cmd = 'abaqus'; end
if nargin < 4 || isempty(simDir),     simDir = 'WorkDir_Abaqus'; end
if nargin < 5,                         res_name = ''; end
if nargin < 6 || isempty(doPlot),      doPlot = false; end

Lx = params.Lx;
Ly = params.Ly;
mesh_size = params.mesh_size;

% ----------------- run 4 corrector jobs -----------------
Jobs  = {'Uniaxialx', 'Uniaxialy', 'Shear12', 'Shear21'};
Loads = [0.01, 0,    0,    0; ...
         0,    0,    0,    0.01; ...
         0,    0.01, 0,    0; ...
         0,    0,    0.01, 0];

Results = struct();

for k = 1:4
    % Create Work Directory
    if ~exist(simDir, 'dir'), mkdir(simDir); end
    AllFunctions.Mat_5A_safeCleanRunDir(simDir);

    job_name = Jobs{k};
    H_val    = Loads(k,:);

    fprintf('\n>> Processing RUC Case %d: %s\n', k, job_name);

    % Build parameters for this RUC run
    parms = struct();
    parms.PBC       = PBC_VALUE;
    parms.H11       = H_val(1);
    parms.H12       = H_val(2);
    parms.H21       = H_val(3);
    parms.H22       = H_val(4);
    parms.NX        = 1;
    parms.NY        = 1;
    parms.JOB       = job_name;
    parms.Lx        = Lx;
    parms.Ly        = Ly;
    parms.mesh_size = mesh_size;

    % Generate python script
    out_py_name  = [job_name, '_run.py'];
    full_py_path = fullfile(simDir, out_py_name);

    fprintf('>> Generating script using AllFunctions.Mat_7A_updatePyFromParams...\n');
    AllFunctions.Mat_7A_updatePyFromParams(parms, pyFile, full_py_path);

    % Run Abaqus (Windows cmd style, same as your script)
    fprintf('>> Running Abaqus Job: %s\n', job_name);
    cmd = sprintf('cmd /c "cd /d ""%s"" && %s cae noGUI=""%s"""', ...
        fullfile(pwd, simDir), abaqus_cmd, out_py_name);
    status = system(cmd);
    if status ~= 0
        error('Abaqus simulation failed. Check .log files in %s', simDir);
    end

    % Load CSV outputs
    csv_path  = fullfile(simDir, sprintf('%s_NodalData.csv', job_name));
    conn_path = fullfile(simDir, sprintf('%s_Connectivity.csv', job_name));

    if ~isfile(csv_path)
        error('RUC simulation failed for %s. Missing file: %s', job_name, csv_path);
    end
    if ~isfile(conn_path)
        error('RUC simulation failed for %s. Missing file: %s', job_name, conn_path);
    end

    data     = readtable(csv_path);
    raw_conn = readmatrix(conn_path, 'NumHeaderLines', 1);
    Connectivity = raw_conn(:, :);

    % Processing (same as your script)
    Y1 = data.X; Y2 = data.Y;
    U1 = data.U_U1; U2 = data.U_U2;
    Node_Labels = data.Node_Label;

    U1_Aff = H_val(1)*Y1 + H_val(2)*Y2;
    U2_Aff = H_val(3)*Y1 + H_val(4)*Y2;

    denom = max(abs(H_val));
    Xi1 = (U1 - U1_Aff) / denom;
    Xi2 = (U2 - U2_Aff) / denom;

    GradChi1 = AllFunctions.Mat_2C_FEM_2D_GaussPoint_gradient_from_NodalPoint( ...
        Y1, Y2, Xi1, Connectivity, Node_Labels);
    GradChi2 = AllFunctions.Mat_2C_FEM_2D_GaussPoint_gradient_from_NodalPoint( ...
        Y1, Y2, Xi2, Connectivity, Node_Labels);

    G11 = GradChi1(:,6);  G12 = GradChi1(:,7);
    G21 = GradChi2(:,6);  G22 = GradChi2(:,7);

    Results(k).Y1 = GradChi1(:,3);
    Results(k).Y2 = GradChi1(:,4);
    Results(k).G11 = G11; Results(k).G12 = G12;
    Results(k).G21 = G21; Results(k).G22 = G22;
    Results(k).NodeLabels = Node_Labels;

    if doPlot
        figure('Color','w','Name',['Fluctuation Gradients: ' job_name], 'Position', [100, 100, 1000, 800]);
        titles = {'\partial \Xi_1 / \partial y_1', '\partial \Xi_1 / \partial y_2', ...
                  '\partial \Xi_2 / \partial y_1', '\partial \Xi_2 / \partial y_2'};
        grads = {G11, G12, G21, G22};
        for sp = 1:4
            subplot(2,2,sp);
            scatter3(Results(k).Y1, Results(k).Y2, grads{sp}, 15, grads{sp}, 'filled');
            title(titles{sp}); colorbar; grid on; view(45,30);
        end
    end
end

% ----------------- assemble localization tensor -----------------
fprintf('\n--- ASSEMBLING LOCALIZATION TENSOR --- \n');

R11 = Results(1); R22 = Results(2); R12 = Results(3); R21 = Results(4);

Y1 = R11.Y1; Y2 = R11.Y2;
nN = numel(Y1);

dv11 = cat(3, [R11.G11, R11.G12], [R11.G21, R11.G22]);
dv22 = cat(3, [R22.G11, R22.G12], [R22.G21, R22.G22]);
dv12 = cat(3, [R12.G11, R12.G12], [R12.G21, R12.G22]);
dv21 = cat(3, [R21.G11, R21.G12], [R21.G21, R21.G22]);

grads_list = {dv11, dv22, dv12, dv21};

% Index mapping for the 4-component basis: [11, 22, 12, 21]
pair  = @(idx) (idx==1)*[1 1] + (idx==2)*[2 2] + (idx==3)*[1 2] + (idx==4)*[2 1];
delta = @(a,b) double(a==b);

Lmat = zeros(nN, 4, 4);

for p = 1:nN
    for mc = 1:4
        dv_node = squeeze(grads_list{mc}(p,:,:))';  % 2x2, dv(i,j)

        for sc = 1:4
            ij = pair(sc); i = ij(1); j = ij(2);

            % I_ij basis (4 independent cases)
            I_ij = (mc==1) * (delta(i,1)*delta(j,1)) ...
                 + (mc==2) * (delta(i,2)*delta(j,2)) ...
                 + (mc==3) * (delta(i,1)*delta(j,2)) ...
                 + (mc==4) * (delta(i,2)*delta(j,1));

            Lmat(p, sc, mc) = I_ij + dv_node(i,j);
        end
    end
end

L = struct();
L.Y1 = Y1; L.Y2 = Y2;
L.NodeLabels = R11.NodeLabels;

% Keep EXACT same field naming/order as your script
L.L11_E11 = Lmat(:,1,1); L.L11_E22 = Lmat(:,1,2); L.L11_E12 = Lmat(:,1,3); L.L11_E21 = Lmat(:,1,4);
L.L12_E11 = Lmat(:,3,1); L.L12_E22 = Lmat(:,3,2); L.L12_E12 = Lmat(:,3,3); L.L12_E21 = Lmat(:,3,4);
L.L21_E11 = Lmat(:,4,1); L.L21_E22 = Lmat(:,4,2); L.L21_E12 = Lmat(:,4,3); L.L21_E21 = Lmat(:,4,4);
L.L22_E11 = Lmat(:,2,1); L.L22_E22 = Lmat(:,2,2); L.L22_E12 = Lmat(:,2,3); L.L22_E21 = Lmat(:,2,4);

% ----------------- optional save -----------------
if ~isempty(res_name)
    save(res_name, 'L');
end
save('LocalizationTensor_RUC.mat','L');

end

