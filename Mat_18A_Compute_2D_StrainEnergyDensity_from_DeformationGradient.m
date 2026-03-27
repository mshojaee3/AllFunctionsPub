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
%J  = F11.*F22 - F12.*F21;
%I1 = 1 + F11.^2 + F12.^2 + F21.^2 + F22.^2;

% ----------------- strain energy -----------------
%W = theta(1).*(J.^(-2/3).*I1 - 3) ...
%  + theta(3).*(J - 1).^2;


%W = ((-3)+((-1).*F12.*F21+F11.*F22).^(-2/3).*(1+F11.^2+F12.^2+F21.^2+ ...
%    F22.^2)).*theta(1)+((-3)+((-1).*F12.*F21+F11.*F22).^(-4/3).*( ...
%    F21.^2+F12.^2.*(1+F21.^2)+(-2).*F11.*F12.*F21.*F22+F22.^2+F11.^2.* ...
%    (1+F22.^2))).*theta(2)+((-1)+(-1).*F12.*F21+F11.*F22).^2.*theta(3) ...
%    +((-3)+((-1).*F12.*F21+F11.*F22).^(-4/3).*(1+(F11.^2+F21.^2).^2+( ...
%    F12.^2+F22.^2).^2)).*theta(4)+((-3)+((-1).*F12.*F21+F11.*F22).^( ...
%    -2).*(1+(F11.^2+F21.^2).^3+(F12.^2+F22.^2).^3)).*theta(5)+((-3)+(( ...
%    -1).*F12.*F21+F11.*F22).^(-8/3).*(1+(F11.^2+F21.^2).^4+(F12.^2+ ...
%    F22.^2).^4)).*theta(6)+((-3)+((-1).*F12.*F21+F11.*F22).^(-4).*(1+( ...
%    F11.^2+F21.^2).^6+(F12.^2+F22.^2).^6)).*theta(7)+((-3)+((-1).* ...
%    F12.*F21+F11.*F22).^(-2/3).*(1+F11.^2+F12.^2+F21.^2+F22.^2)).^2.* ...
%    theta(8)+((-3)+((-1).*F12.*F21+F11.*F22).^(-2/3).*(1+F11.^2+ ...
%    F12.^2+F21.^2+F22.^2)).*((-3)+((-1).*F12.*F21+F11.*F22).^(-4/3).*( ...
%    F21.^2+F12.^2.*(1+F21.^2)+(-2).*F11.*F12.*F21.*F22+F22.^2+F11.^2.* ...
%    (1+F22.^2))).*theta(9)+((-3)+((-1).*F12.*F21+F11.*F22).^(-4/3).*( ...
%    F21.^2+F12.^2.*(1+F21.^2)+(-2).*F11.*F12.*F21.*F22+F22.^2+F11.^2.* ...
%    (1+F22.^2))).^2.*theta(10)+((-1)+(-1).*F12.*F21+F11.*F22).^4.* ...
%    theta(11);


J   = F11.*F22 - F12.*F21;

C11 = F11.^2 + F21.^2;
C22 = F12.^2 + F22.^2;
C12 = F11.*F12 + F21.*F22;

% plane strain invariants (3D embedding with F33 = 1)
I1  = 1 + C11 + C22;
I2  = C11 + C22 + C11.*C22 - C12.^2;   % same as C11 + C22 + J.^2

I1b = J.^(-2/3) .* I1;
I2b = J.^(-4/3) .* I2;

% cubic / anisotropic invariant-like terms
J7b = J.^(-4/3) .* (1 + C11.^2 + C22.^2);
J8b = J.^(-2)   .* (1 + C11.^3 + C22.^3);
J9b = J.^(-8/3) .* (1 + C11.^4 + C22.^4);

W = theta(1).*(I1b - 3) ...
  + theta(2).*(I2b - 3) ...
  + theta(3).*(I1b - 3).^2 ...
  + theta(4).*(I1b - 3).*(I2b - 3) ...
  + theta(5).*(I2b - 3).^2 ...
  + theta(6).*(J - 1).^2 ...
  + theta(7).*(J - 1).^4 ...
  + theta(8).*(J7b - 3) ...
  + theta(9).*(J8b - 3) ...
  + theta(10).*(J9b - 3);


end
