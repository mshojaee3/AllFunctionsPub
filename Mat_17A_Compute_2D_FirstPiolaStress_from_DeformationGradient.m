function [P11,P12,P21,P22] = Mat_17A_Compute_2D_FirstPiolaStress_from_DeformationGradient(params)
% firstPiola
% Computes 2D First Piola-Kirchhoff stress components from
% deformation gradient components and material parameters.
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
%   P11, P12, P21, P22 : First Piola stress components

% ----------------- validate -----------------
assert(isstruct(params), 'params must be a struct.');
assert(isfield(params,'theta'), 'params.theta must be provided.');
assert(isfield(params,'F11'),   'params.F11 must be provided.');
assert(isfield(params,'F12'),   'params.F12 must be provided.');
assert(isfield(params,'F21'),   'params.F21 must be provided.');
assert(isfield(params,'F22'),   'params.F22 must be provided.');

theta = params.theta;
F11   = params.F11;
F12   = params.F12;
F21   = params.F21;
F22   = params.F22;

J = F11*F22 - F12*F21;

P11 = theta(1)*F11 + theta(3)*J*F22;
P12 = theta(1)*F12 - theta(3)*J*F21;
P21 = theta(1)*F21 - theta(3)*J*F12;
P22 = theta(1)*F22 + theta(3)*J*F11;

end
