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

J  = F11.*F22 - F12.*F21;
I1 = 1 + F11.^2 + F12.^2 + F21.^2 + F22.^2;

P11 = (2.*F11.*J.^(-2/3) - (2/3).*F22.*J.^(-5/3).*I1).*theta(1) ...
    + 2.*F22.*(J-1).*theta(3);

P12 = (2.*F12.*J.^(-2/3) + (2/3).*F21.*J.^(-5/3).*I1).*theta(1) ...
    - 2.*F21.*(J-1).*theta(3);

P21 = (2.*F21.*J.^(-2/3) + (2/3).*F12.*J.^(-5/3).*I1).*theta(1) ...
    - 2.*F12.*(J-1).*theta(3);

P22 = (2.*F22.*J.^(-2/3) - (2/3).*F11.*J.^(-5/3).*I1).*theta(1) ...
    + 2.*F11.*(J-1).*theta(3);

end
