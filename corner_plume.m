function yDot = corner_plume(z,X,const,options,ambient)
%QUARTER_PLUME
%
% This function is called by BPTmodel. It is not intended to be used as a stand-alone function.
%
% See also BPTmodel

% INPUTS
%   - z = depth
%   - X(1) = the plume radius (m)
%   - X(2) = the plume velocity (m s^-1)
%   - X(3) = the plume temperature (C)
%   - X(4) = the plume salinity (psu)
%
%   - const = structure containing preset constants in the following fields: 
%       - g = gravity accel
%       - rho0 = reference density
%       - ci / cw = heat capacity of ice / water
%       - lamba* = coefficients for linearized freezing temp equation
%       - L = latent heat of fusion
%   - options = structure containing user-defined (or default) constants in the following fields:
%       - alpha = entrainment coefficient
%       - Cd = drag coefficient
%       - gammaT = turbulent transfer coefficient for temp
%       - gammaS = turbulent transfer coefficient for salt
%       - iceT = temperature of ice
%   - ambient = structure containing user-defined values for ambient conditions
%       - temp / salt / depth = properties of ambient water
%       - u_horiz = horizontal velocity that enhances melt rates

% OUTPUTS
%   - yDot(1) = d/dz of plume radius (m)
%   - yDot(2) = d/dz of plume velocity (m s^-1)
%   - yDot(3) = d/dz of plume temperature (C)
%   - yDot(4) = d/dz of plume salinity (psu)

% FUNCTIONS called:
%       melt_calc(...)  
    
% DESCRIPTION:
%       this code solves BPT assuming a QUARTER CONICAL plume geometry with two faces in contact with ice
yDot = zeros(4,1);

% calculate melt rate with calc_melt()
wmelt = sqrt(X(2).^2 + ambient.u_horiz.^2);

[melt, Tb, Sb] = melt_calc(z, wmelt, X(3), X(4), lambda1=const.lambda1, lambda2=const.lambda2, lambda3=const.lambda3, gammaT=options.gammaT, gammaS=options.gammaS, Cd=options.Cd, iceT=options.iceT, L=const.L, ci=const.ci, cw=const.cw);

%ODEs for radius, velocity, temperature, and salinity. 
yDot(1)=2*options.alpha+8*melt/(pi*X(2))-X(1)*const.g*(interp1(ambient.depth,ambient.density,z)-sw_pden(X(4),X(3),abs(z),0))/(2*X(2)*X(2)*const.rho0)+4*options.Cd/pi;

yDot(2)=-2*options.alpha*X(2)/X(1)-8*melt/(pi*X(1))+const.g*(interp1(ambient.depth,ambient.density,z)-sw_pden(X(4),X(3),abs(z),0))/(X(2)*const.rho0)-8*options.Cd*X(2)/(pi*X(1));

yDot(3)=2*options.alpha*(interp1(ambient.depth,ambient.temp,z)-X(3))/X(1)+8*melt*(Tb-X(3))/(pi*X(1)*X(2))-8*options.gammaT*sqrt(options.Cd)*(X(3)-Tb)/(pi*X(1))*wmelt/X(2);

yDot(4)=2*options.alpha*(interp1(ambient.depth,ambient.salt,z)-X(4))/X(1)+8*melt*(Sb-X(4))/(pi*X(1)*X(2))-8*options.gammaS*sqrt(options.Cd)*(X(4)-Sb)/(pi*X(1))*wmelt/X(2);

end
