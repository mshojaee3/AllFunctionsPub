function W = Mat_18A_Compute_2D_StrainEnergyDensity_from_DeformationGradient(params)
% strainEnergy
% Computes strain energy density from deformation gradient components
% and material parameters stored in one input struct.
%
% INPUT
%   params : struct with fields
%            params.theta
%            params.F11
%            params.F12
%            params.F21
%            params.F22
%
% OUTPUT
%   W : strain energy density

% ----------------- validate -----------------
assert(isstruct(params), 'params must be a struct.');
assert(isfield(params,'theta'), 'params.theta must be provided.');
assert(isfield(params,'F11'), 'params.F11 must be provided.');
assert(isfield(params,'F12'), 'params.F12 must be provided.');
assert(isfield(params,'F21'), 'params.F21 must be provided.');
assert(isfield(params,'F22'), 'params.F22 must be provided.');

theta = params.theta;
F11   = params.F11;
F12   = params.F12;
F21   = params.F21;
F22   = params.F22;

% ----------------- invariants -----------------
J  = F11.*F22 - F12.*F21;
I1 = 1 + F11.^2 + F12.^2 + F21.^2 + F22.^2;

% ----------------- strain energy -----------------
W = theta(1).*(J.^(-2/3).*I1 - 3) ...
  + theta(3).*(J - 1).^2;

end
