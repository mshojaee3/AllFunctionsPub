function Mat_11A_plot_2D_H_abaqus_vs_predicted(parms)
% Plot-only. Uses parms.plotdata.(coor, H_abaqus, H_pred)

% ---- check required packed data ----

pd = parms.plotdata;

req = {'coor','H_abaqus','H_pred'};
for i = 1:numel(req)
    if ~isfield(pd, req{i})
        error('parms.plotdata.%s is missing.', req{i});
    end
end

coor     = pd.coor;
H_abaqus = pd.H_abaqus;
H_pred   = pd.H_pred;

% ---- optional settings ----
tag = '';
if isfield(pd,'nameTag') && ~isempty(pd.nameTag)
    tag = char(pd.nameTag);
end


x_gp = coor(:,1);
y_gp = coor(:,2);

% ---- Figure 1: Abaqus ----
namesA = {'H11^{Abaqus}','H12^{Abaqus}','H21^{Abaqus}','H22^{Abaqus}'};
figure('Color','w','Name',['H_Abaqus ' tag],'Position',[100 100 1100 850]);
for k = 1:4
    subplot(2,2,k);
    scatter3(x_gp, y_gp, H_abaqus{k}, 15, H_abaqus{k}, 'filled');
    title(namesA{k}); xlabel('X'); ylabel('Y'); zlabel(namesA{k});
    colorbar; grid on; view(45,30);
end

% ---- Figure 2: Predicted ----
namesP = {'H11^{Pred}','H12^{Pred}','H21^{Pred}','H22^{Pred}'};
figure('Color','w','Name',['H_Predicted ' tag],'Position',[150 130 1100 850]);
for k = 1:4
    subplot(2,2,k);
    scatter3(x_gp, y_gp, H_pred{k}, 15, H_pred{k}, 'filled');
    title(namesP{k}); xlabel('X'); ylabel('Y'); zlabel(namesP{k});
    colorbar; grid on; view(45,30);
end
end
