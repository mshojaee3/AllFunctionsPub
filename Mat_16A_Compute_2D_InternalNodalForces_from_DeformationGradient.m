function InternalNodalForce = Mat_16A_Compute_2D_InternalNodalForces_from_DeformationGradient(params)
% ComputeInternalNodalForces
% Computes internal nodal force vector using deformation-gradient maps,
% IGA basis data, and material parameters stored in one params struct.
%
% INPUT
%   params : struct with fields
%            params.theta
%            params.F11_map
%            params.F12_map
%            params.F21_map
%            params.F22_map
%            params.IsoGeo
%            params.Length_x
%            params.Length_y
%
% OUTPUT
%   InternalNodalForce : internal nodal force vector

theta    = params.theta;
F11_map  = params.F11_map;
F12_map  = params.F12_map;
F21_map  = params.F21_map;
F22_map  = params.F22_map;
IsoGeo   = params.IsoGeo;
Length_x = params.Length_x;
Length_y = params.Length_y;

x_gp = IsoGeo.x_gp;
y_gp = IsoGeo.x_gp;

w_x = IsoGeo.w_gp;
w_y = IsoGeo.w_gp;

N_gp   = IsoGeo.N;
M_gp   = IsoGeo.N;
d1N_gp = IsoGeo.dN;
d1M_gp = IsoGeo.dN;

Nu_x = IsoGeo.Nu_x;
Nu_y = IsoGeo.Nu_x;
n_gp = IsoGeo.n_gp;
n_xi = IsoGeo.n_xi;
n_eta = IsoGeo.n_xi;

xi_node  = linspace(0,1,Nu_x+1);
eta_node = linspace(0,1,Nu_y+1);

num_total_gp_x = Nu_x * n_gp;
num_total_gp_y = Nu_y * n_gp;

gp_P11 = zeros(num_total_gp_y, num_total_gp_x);
gp_P12 = zeros(num_total_gp_y, num_total_gp_x);
gp_P21 = zeros(num_total_gp_y, num_total_gp_x);
gp_P22 = zeros(num_total_gp_y, num_total_gp_x);

Fx = zeros(n_xi, n_eta);
Fy = zeros(n_xi, n_eta);

for ex = 1:Nu_x
    for ey = 1:Nu_y
        for gx = 1:n_gp
            for gy = 1:n_gp

                ix = (ex-1)*n_gp + gx;
                iy = (ey-1)*n_gp + gy;

                X = x_gp(ix) * Length_x;
                Y = y_gp(iy) * Length_y;

                F11 = F11_map(X,Y);
                F12 = F12_map(X,Y);
                F21 = F21_map(X,Y);
                F22 = F22_map(X,Y);

                wxy = w_x(ix) * w_y(iy) * Length_x * Length_y * ...
                      (xi_node(ex+1)-xi_node(ex)) * ...
                      (eta_node(ey+1)-eta_node(ey)) / 4;

                paramsPiola = struct();
                paramsPiola.theta = theta;
                paramsPiola.F11   = F11;
                paramsPiola.F12   = F12;
                paramsPiola.F21   = F21;
                paramsPiola.F22   = F22;

                [P11,P12,P21,P22] = firstPiola(paramsPiola);

                gp_P11(iy,ix) = P11;
                gp_P12(iy,ix) = P12;
                gp_P21(iy,ix) = P21;
                gp_P22(iy,ix) = P22;

                for a = 1:n_xi
                    for b = 1:n_eta
                        dRdx = (1/Length_x) * d1N_gp(a,ix) * M_gp(b,iy);
                        dRdy = (1/Length_y) * N_gp(a,ix)   * d1M_gp(b,iy);

                        Fx(a,b) = Fx(a,b) + (P11*dRdx + P12*dRdy) * wxy;
                        Fy(a,b) = Fy(a,b) + (P21*dRdx + P22*dRdy) * wxy;
                    end
                end
            end
        end
    end
end

Fx_int = Fx(2:end-1, 2:end-1);
Fy_int = Fy(2:end-1, 2:end-1);

InternalNodalForce = [Fx_int(:); Fy_int(:)];

end
