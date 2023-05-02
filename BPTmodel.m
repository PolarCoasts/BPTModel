function experiment = BPTmodel(aD,aT,aS,GL,Q,options)

arguments
   aD (1,:) {mustBeLessThanOrEqual(aD,0)}                               % ambient profile depths (m)
   aT (1,:) {mustBeEqualSize(aD,aT)}                                    % ambient profile of temperature (degrees C)
   aS (1,:) {mustBeEqualSize(aD,aS)}                                    % ambient profile of salinity (psu)
   GL (1,1) {mustBeLessThan(GL,0)}                                      % grounding line depth (m)
   Q (1,1) {mustBeGreaterThan(Q,0)}                                     % subglacial discharge flux (m^3/s)
   options.iceT (1,1) double=-10                                        % ice temperature (degrees C)
   options.W (1,1) double=100                                           % outlet width (m)
   options.uh (1,1) double=0                                            % horizontal velocity relevant for melting (m/s)
   options.ui {mustBeValueOrChar(options.ui,"balance")}='balance'       % initial velocity (m/s), 'balance' of momentum and buoyancy or set value
   options.alpha (1,1) double=0.1                                       % entrainment coefficient
   options.Cd (1,1) double=2.5e-3                                       % drag coefficient
   options.gammaT (1,1) double=0.022                                    % thermal transfer coefficient
   options.gammaS (1,1) double=0.00062                                  % haline transfer coefficient
   options.intMethod char {mustBeMember(options.intMethod,{'linear','nearest','next','previous','pchip','cubic','v5cubic','makima','spline'})}='linear'  % method for interpolating ambient profile
   options.extValue_T {mustBeValueOrChar(options.extValue_T,"extrap")}='extrap' % to extend ambient profile to full plume range, extrapolate or set to specific value
   options.extValue_S {mustBeValueOrChar(options.extValue_S,"extrap")}='extrap' % to extend ambient profile to full plume range, extrapolate or set to specific value
end
% INPUTS:
%     required:
%     - aD / aT / aS = depth, temperature, and salinity of ambient water
%     - GL = grounding line depth
%     - Q = subglacial discharge flux
%     
%     optional (default):
%     - iceT (-10) = temperature of ice
%     - W (100) = discharge outlet width
%     - uh (0) = horizontal velocity to be used for calculating melt
%     - ui ('balance') = initial upwelling velocity (specify value or 'balance')
%     - alpha (0.1) = entrainment coefficient
%     - Cd (2.5e-3) = drag coefficient
%     - gammaT (0.022) = thermal transfer coefficient
%     - gammaS (0.00062) = haline transfer coefficient
%     - intMethod ('linear') = interpolation method for ambient profile (methods of interpolation for interp1)
%     - extValue_T ('extrap') = value to fill ambient temperature profile with to extend to full water depth (or specify 'extrap')
%     - extValue_S ('extrap') = value to fill ambient salinity profile with to extend to full water depth (or specify 'extrap')
% 
% examples: experiment = BPTline(aD,aT,aS,GL,Q)  *minimum inputs
%           experiment = BPTline(aD,aT,aS,GL,Q,W=200,alpha=0.11) *setting two optional values
% 
% OUTPUT (structure with four fields):
%     plume:
%     - depth / temp / salt / density = properties of upwelling plume
%     - w = vertical velocity of plume
%     - melt = melt rate
%     - neutralDensity = depth where plume density equals ambient density
%     - maximumH = depth of maximum height reached by plume (momentum=0 or surface)
%     - radius = thickness of plume (perpendicular to ice)
%     - maximumMD = depth of maximum melt rate
%     - volumeFlux = volume flux of plume
%     - momentumFlux = momentum flux of plume (**units need to be checked**)
%     - units = list of units for plume variables
% 
%     ambient:
%     - depth / temp / salt = input profiles interpolated to 0.1 m increments and extended to full water column
%     - u_horiz = input horizontal velocity
%     - density = density of ambient profile
%     - N2 = buoyancy frequency of ambient profile
%     - units = list of units for ambient variables
%        
%     glacier:
%     - GL = grounding line depth
%     - outletW = discharge outlet width
%     - Q = subglacial discharge flux
%     - iceTemp = temperature of ice
%     - units = list of units for glacier variables
%     
%     coeff:
%     - alpha = entrainment coefficient
%     - Cd = drag coefficient
%     - gammaT = thermal transfer coefficient
%     - gammaS = haline transfer coefficient
% 
% FUNCTIONS called:
%       line_plume(...)
%       melt_calc(...)
   
% calculate discharge as m2/s (total discharge in m3 / width of plume)
qsg = Q/options.W;

% interpolate ambient profile to 0.1 m depth increments, extrapolate/fill to full depth range
ii=find(~isnan(aT) & ~isnan(aS));
aD_good=aD(ii); aT_good=aT(ii); aS_good=aS(ii);
depth=0:-.1:GL;
ambientTemp=interp1(aD_good,aT_good,depth,options.intMethod,options.extValue_T);
ambientSalt=interp1(aD_good,aS_good,depth,options.intMethod,options.extValue_S);

%calculate ambient density profile 
ambientRho = sw_pden(ambientSalt,ambientTemp,abs(depth),0);

% Include interpolated input variables in output structure
experiment.ambient.depth = depth;
experiment.ambient.temp = ambientTemp;
experiment.ambient.salt = ambientSalt;
experiment.ambient.density = ambientRho;
experiment.ambient.u_horiz = options.uh;

%% Define Constants

const.g = 9.81; %gravitational acceleration (m s^-2)

const.cw = 3974; %heat capacity of seawater (J kg^-1 C^-1)
const.ci = 2009.0; %heat capacity of ice  (J kg^-1 C^-1)
const.L = 335000; %latent heat of fusion (J kg^-1)

const.rho0 = 1028; %reference density (kg m^-3) 

const.lambda1 = -0.0573; %variation of freezing point with salinity (C psu^-1)
const.lambda2 = 0.0832; %freezing point offset (C)
const.lambda3 = 0.000761; %variation of freezing point with depth (C m^-1)

%% Initial conditions of plume at grounding line

ti = sw_fp(0,abs(GL)); %initial plume starts at freezing point
si = 0.0001; % initial plume salinity  (note: needs > 0 for integration to converge)

surface = 0;
nz = length(ambientTemp);

% find ambient rho and plume rho at grounding line --> calculate g' (gp)
iGL = find(abs(depth)==abs(GL),1);
ambientRho_GL = ambientRho(iGL);
plumeRho_GL = sw_pden(si,ti,abs(GL),0);
gp = (ambientRho_GL-plumeRho_GL)*const.g/const.rho0;

% initial velocity of plume
if ismember(options.ui,"balance")
    ui = (gp*qsg/options.alpha)^(1/3); %balance of momentum and buoyancy for a line
else
    ui=options.ui;
end

bi = qsg / ui; %plume radius or thickness in cross-terminus direction (m)

%% Solve coupled system of ODEs 
% using a 4th order MATLAB integrator
%    z = depth
%    X is m x 4 with columns of plume: width, velocity, T , S

[z,X]=ode45(@(Z,x) line_plume(Z,x,const,options,experiment.ambient),[GL,surface],[bi,ui,ti,si]);

%% Compute additional variables from solution of ODEs

[m n] = size(X);
A = zeros(m,4);

% quantities from solving ODEs
    % X(:,1) = radius / width of plume (m)
    % X(:,2) = velocity of plume (m/s)
    % X(:,3) = plume temperature (C)
    % X(:,4) = plume salinity (psu)

% quantities to calculate in the loop here
    % A(:,1) = volume flux (m2/s)
    % A(:,2) = momentum flux (m4/s)
    % A(:,3) = melt rate (m/day)
    % A(:,4) = plume density (kg/m3)

for i=1:m   
    %volume flux (m2/s)
    A(i,1)=X(i,1)*X(i,2);
    
    % calculate total velocity that drives melting: upwelling + horizontal 
    uTot = sqrt(X(i,2).^2 + options.uh.^2); 

    %melt rate in m/s --> m/day
    [melt_ms,~,~] = melt_calc(z(i), uTot, X(i,3), X(i,4), lambda1=const.lambda1, lambda2=const.lambda2, lambda3=const.lambda3, gammas=options.gammaS, gammaT=options.gammaT, Cd=options.Cd, iceT=options.iceT, L=const.L, ci=const.ci, cw=const.cw);
    A(i,3)=24*60*60*melt_ms;
    
    %plume density (kg m^-3)
    A(i,4)=sw_pden(X(i,4),X(i,3),abs(z(i)),0);
    
    %momentum flux (m^4 s^-1 - not right units?) - weird, this has a rho here too??
    A(i,2)=X(i,1)*X(i,2)*X(i,2)*A(i,4);   
end

%% Interpolate variables onto specified vertical grid and give final names
experiment.plume.depth = depth;
experiment.plume.radius = interp1(z,X(:,1),depth);
experiment.plume.w = interp1(z,X(:,2),depth); 
experiment.plume.temp = interp1(z,X(:,3),depth);
experiment.plume.salt = interp1(z,X(:,4),depth); 
experiment.plume.density = interp1(z,A(:,4),depth);
experiment.plume.volumeFlux = interp1(z,A(:,1),depth);
experiment.plume.momentumFlux = interp1(z,A(:,2),depth); 
experiment.plume.melt = interp1(z,A(:,3),depth);

%find depth of maximum melt
mi = find(experiment.plume.melt == max(experiment.plume.melt),1);
experiment.plume.maximumMD = depth(mi); 

% find terminal level based on density
tl = find(experiment.plume.density <= ambientRho,1);
experiment.plume.neutralDensity = depth(tl);

%find where integration stops
wIndex = find(isnan(experiment.plume.w));

%if integration did not stop, then plume reached surface
if isempty(wIndex)
    wIndex = 1;
end

%find maximum height (which is above terminal depth)
experiment.plume.maximumH = depth(wIndex(end));

%add units for all plume variables
experiment.plume.units={'depth (m)'; 'radius (m)'; 'w (m/s)'; 'temp (degrees C)'; 'salt (psu)'; 'density (kg/m^3)'; 'volumeFlux (m^3/s)'; 'momentumFlux (need to check)'; 'melt (m/day)'; 'maximumMD (m)'; 'neutralDensity (m)'; 'maximumH (m)'};

%% calculate N2 of ambient profile
experiment.ambient.N2 = real(sqrt(const.g/const.rho0 * diff(ambientRho))); 

%add units for all ambient variables
experiment.ambient.units={'depth (m)'; 'temp (degrees C)'; 'salt (psu)'; 'density (kg/m^3)'; 'u_horiz (m/s)'; 'N2 (s^-1)'};

%% Return glacier inputs
experiment.glacier.GL = GL;
experiment.glacier.outletW = options.W;
experiment.glacier.Q = Q;
experiment.glacier.iceTemp = options.iceT;
experiment.glacier.units={'GL (m)'; 'outletW (m)'; 'Q (m^3/s)'; 'iceTemp (degrees C)'};

%% Return coefficients used
experiment.coeff.alpha = options.alpha;
experiment.coeff.Cd = options.Cd;
experiment.coeff.gammaT = options.gammaT;
experiment.coeff.gammaS = options.gammaS;

end
%% Validation functions
function mustBeEqualSize(a,b)
    if ~isequal(size(a),size(b))
        eid='Size:notEqual';
        msg='Vectors for ambient properties must all have the same dimensions';
        throwAsCaller(MException(eid,msg))
    end
end

function mustBeValueOrChar(a,string)
    if ~ismember(a,string) && ~isscalar(a)
        eid='Value:Error';
        msg="extValue must be 'extrapolate' or a constant value";
        throwAsCaller(MException(eid,msg))
    end
end
    