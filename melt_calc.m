function [melt, Tb, Sb] = melt_calc(z, u, T, S, const, options)

arguments
    z (1,:) {mustBeLessThanOrEqual(z,0)}        % depth
    u (1,:) {mustBeEqualSize(u,z)}              % velocity that drives melting
    T (1,:) {mustBeEqualSize(T,z)}              % temperature outside the BL (e.g. in plume)
    S (1,:) {mustBeEqualSize(S,z)}              % salinity outside the BL (e.g. in plume)
    const.cw (1,1) double=3974                  % heat capacity of seawater (J kg^-1 C^-1)
    const.ci (1,1) double=2009                  % heat capacity of ice (J kg^-1 C^-1)
    const.L (1,1) double=335000;                % latent heat of fusion (J kg^-1)
    const.lambda1 (1,1) double=-0.0573          % variation of freezing point with salinity (C psu^-1)
    const.lambda2 (1,1) double=0.0832           % freezing point offset (C)
    const.lambda3 (1,1) double=0.000761         % variation of freezing point with depth (C m^-1)
    options.iceT (1,1) double=-10               % ice temperature (degrees C)
    options.Cd (1,1) double=2.5e-3              % drag coefficient
    options.gammaT (1,1) double=0.022           % thermal transfer coefficient
    options.gammaS (1,1) double=0.00062         % haline transfer coefficient
end

% INPUTS
%   required:
%   - z = depth
%   - u = velocity that drives melting
%   - T = temperature outside the BL (e.g. in plume)
%   - S = temperature outside the BL (e.g. in plume)
%   
%   constants (default): 
%   - ci / cw = heat capacity of ice / water
%   - L = latent heat of fusion
%   - lamba* = coefficients for linearized freezing temp equation
% 
%   options (default):
%   - gammaS = turbulent transfer coefficient for salt 
%   - gammaT = turbulent transfer coefficient for temp
%   - Cd = drag coefficient
%   - iceT = temperature of ice
%
% OUTPUTS
%   - melt = melt rate in m/s
%   - Tb = temperature at ice-ocean boundary
%   - Sb = salinity at ice-ocean boundary
%
% examples: [melt, Tb, Sb] = melt_calc(z, u, T, S)  *minimum inputs
%           [melt, Tb, Sb] = melt_calc(z, u, T, S,iceT=0,gammaT=0.044) *setting two optional values
%
% DESCRIPTION:
%       3-Equation Melt Parameterization (e.g. Jenkins 2011) 
%       solve quadratic equation for Sb, Tb, and melt
%       based on z, u, T, S and coefficients


% define a, b, c to solve quadratic equation for Sb
a = const.lambda1 * (options.gammaT * const.cw -options.gammaS*const.ci);

b = options.gammaS*const.ci*(const.lambda1 * S - const.lambda2-const.lambda3*z + ...
    options.iceT - (const.L/const.ci)) - options.gammaT*const.cw*(T-const.lambda2-const.lambda3*z);

c = options.gammaS * S * (const.ci * (const.lambda2 + const.lambda3*z-options.iceT)+const.L);

% calculate Sb from a, b, c
Sb = (1/(2*a)) * (-b-sqrt(b^2 - 4*a*c));
% calculate Tb from linearized freezing temp. equation
Tb = const.lambda1*Sb+const.lambda2+const.lambda3*z;

% calculate melt
melt = options.gammaS * sqrt(options.Cd) * u * (S-Sb)/Sb;

end

%% Validation functions
function mustBeEqualSize(a,b)
    if ~isequal(size(a),size(b))
        eid='Size:notEqual';
        msg='Vectors for ambient properties must all have the same dimensions';
        throwAsCaller(MException(eid,msg))
    end
end
