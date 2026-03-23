function Mat_13A_3D_Surface_Plotter(parms)
% General single 3D scatter plot
% Uses parms.plotdata.(x, y, z)
%
% Required:
%   parms.plotdata.x
%   parms.plotdata.y
%   parms.plotdata.z
%
% Optional:
%   parms.plotdata.nameTag
%   parms.plotdata.plotTitle
%   parms.plotdata.zLabel

% ---- check required packed data ----
pd = parms.plotdata;

req = {'x','y','z'};
for i = 1:numel(req)
    if ~isfield(pd, req{i})
        error('parms.plotdata.%s is missing.', req{i});
    end
end

x = pd.x;
y = pd.y;
z = pd.z;

% ---- basic size check ----
if numel(x) ~= numel(y) || numel(x) ~= numel(z)
    error('x, y, and z must have the same number of elements.');
end

% make sure they are column vectors
x = x(:);
y = y(:);
z = z(:);

% ---- optional settings ----
tag = '';
if isfield(pd,'nameTag') && ~isempty(pd.nameTag)
    tag = char(pd.nameTag);
end

plotTitle = '3D Scatter Plot';
if isfield(pd,'plotTitle') && ~isempty(pd.plotTitle)
    plotTitle = char(pd.plotTitle);
end

zLab = 'Z';
if isfield(pd,'zLabel') && ~isempty(pd.zLabel)
    zLab = char(pd.zLabel);
end

% ---- create figure ----
figure('Color','w', ...
       'Name',[plotTitle ' ' tag], ...
       'Position',[100 100 900 700]);

scatter3(x, y, z, 15, z, 'filled');
title([plotTitle ' ' tag]);
xlabel('X');
ylabel('Y');
zlabel(zLab);
colorbar;
grid on;
view(45,30);

end
