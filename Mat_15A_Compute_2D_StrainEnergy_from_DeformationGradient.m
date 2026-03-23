function W_total = Mat_15A_Compute_2D_StrainEnergy_from_DeformationGradient(params)
% Computes total strain energy by integrating the strain energy density
% over the 2D domain using deformation-gradient maps and IGA data.
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
%   W_total : total strain energy

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

Nu_x = IsoGeo.Nu_x;
Nu_y = IsoGeo.Nu_x;
n_gp = IsoGeo.n_gp;

xi_node  = linspace(0,1,Nu_x+1);
eta_node = linspace(0,1,Nu_y+1);

W_total = 0;

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

                energyParams = struct();
                energyParams.theta = theta;
                energyParams.F11   = F11;
                energyParams.F12   = F12;
                energyParams.F21   = F21;
                energyParams.F22   = F22;

                W = AllFunctions.Mat_18A_Compute_2D_StrainEnergyDensity_from_DeformationGradient(energyParams);
                W_total = W_total + W * wxy;

            end
        end
    end
end

end
