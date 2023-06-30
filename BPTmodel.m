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
   options.type char {mustBeMember(options.type,{'line','point','ambient'})}='line' % run model as line- or point-plume, ambient melt runs as stacked line-plumes
end
% INPUTS:
%     required:
%     - aD / aT / aS = depth, temperature, and salinity of ambient water
%     - GL = grounding line depth
%     - Q = subglacial discharge flux (type='ambient' requires small discharge to initiate plume, ex: 1e-10 with W=1)
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
% examples: experiment = BPTmodel(aD,aT,aS,GL,Q)  *minimum inputs
%           experiment = BPTmodel(aD,aT,aS,GL,Q,type='point',Cd=2e-3,alpha=0.11) *setting two optional values, run as point plume
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
%     - momentumFlux = momentum flux of plume
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
%       line_plume(...) or point_plume(...)
%       melt_calc(...)
if strcmp(options.type,'line') || strcmp(options.type,'ambient')
    % calculate discharge as m2/s (discharge per unit outlet width)
    qsg = Q/options.W;
elseif strcmp(options.type,'point')
    % keep discharge as m3/s
    qsg = Q;
end

% interpolate ambient profile to 0.1 m depth increments, extrapolate/fill to full depth range
ii=find(~isnan(aT) & ~isnan(aS));
aD_good=aD(ii); aT_good=aT(ii); aS_good=aS(ii);
depth=round(0:-.1:GL,1);
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
if strcmp(options.type,'line') || strcmp(options.type,'ambient')
    if ismember(options.ui,"balance")
        ui = (gp*qsg/options.alpha)^(1/3); %balance of momentum and buoyancy for a line
    else
        ui=options.ui;
    end
    bi = qsg / ui; %plume radius or thickness in cross-terminus direction (m)
elseif strcmp(options.type,'point')
    if ismember(options.ui,"balance")
        ui = 2/pi*(5*pi^2*gp/(32*options.alpha))^(2/5)*qsg^(1/5);
    else
        ui=options.ui;
    end
    bi = sqrt(2*qsg / (pi*ui)); %plume radius (m)
end
%% Solve coupled system of ODEs 
% using a 4th order MATLAB integrator
%    z = depth
%    X is m x 4 with columns of plume: width, velocity, T , S

if strcmp(options.type,'line')
    [z,X]=ode45(@(Z,x) line_plume(Z,x,const,options,experiment.ambient),[GL,surface],[bi,ui,ti,si]);
elseif strcmp(options.type,'point')
    [z,X]=ode45(@(Z,x) point_plume(Z,x,const,options,experiment.ambient),[GL,surface],[bi,ui,ti,si]);
elseif strcmp(options.type,'ambient')
    di=GL; % initial depth of bottom plume is grounding line depth
    count=0; maxH=[]; nD=[];
    while di<0
        count=count+1;
        [tz,tX]=ode45(@(Z,x) line_plume(Z,x,const,options,experiment.ambient),[di,surface],[bi,ui,ti,si]);
        az(count)={tz};
        aX(count)={tX};
        %recalculate initial conditions for next plume
        di=round(tz(end)+.1,1);
        ti=sw_fp(0,abs(di));
        gp=(ambientRho(abs(depth)==abs(di))-sw_pden(si,ti,abs(di),0))*const.g/const.rho0;
        ui=(gp*qsg/options.alpha)^(1/3);
        bi=qsg/ui;
        % store max plume height for each plume
        maxH=[maxH di-.1];
        % return some progress status to user
        formatSpec='Plume %2.0f outflows at %4.1f m\n  ';
        fprintf(formatSpec,count,di-.1)
        if di<0
            formatSpec2='    Solving for next plume ... \n';
        else
            formatSpec2='    Plume integration complete. Preparing final output...\n';
        end
        fprintf(formatSpec2)
    end
end

%% Interpolate to specified grid 
if strcmp(options.type,'ambient')
    nD=[];
    for i=1:count
        ar(i,:)=interp1(az{i},aX{i}(:,1),depth);
        aw(i,:)=interp1(az{i},aX{i}(:,2),depth);
        at(i,:)=interp1(az{i},aX{i}(:,3),depth);
        as(i,:)=interp1(az{i},aX{i}(:,4),depth);
        % find neutral density for each plume
        arho=sw_pden(as(i,:),at(i,:),depth,0);
        nDi=find(arho<=ambientRho,1);
        nD=[nD depth(nDi)];
    end
    % combine ambient melt plumes into one profile
    radius=sum(ar,'omitnan');
    w=sum(aw,'omitnan');
    temp=sum(at,'omitnan');
    salt=sum(as,'omitnan');
else
    radius=interp1(z,X(:,1),depth);
    w=interp1(z,X(:,2),depth);
    temp=interp1(z,X(:,3),depth);
    salt=interp1(z,X(:,4),depth);
end

%% Compute additional variables from solution of ODEs

% define a geometric factor to account for differing plume areas
if strcmp(options.type,'line') || strcmp(options.type,'ambient')
    geom = options.W;
elseif strcmp(options.type,'point')
    geom = pi.*radius./2;
end

% volume flux (m3/s)
Fv=geom.*radius/2;

% calculate total velocity that drives melting: upwelling + horizontal
uTot=sqrt(w.^2 + options.uh.^2);

% melt rate (m/s --> m/day)
[melt_ms,~,~] = melt_calc(depth, uTot, temp, salt, lambda1=const.lambda1, lambda2=const.lambda2, lambda3=const.lambda3, gammaS=options.gammaS, gammaT=options.gammaT, Cd=options.Cd, iceT=options.iceT, L=const.L, ci=const.ci, cw=const.cw);
melt=24*60*60*melt_ms;

%plume density (kg m^-3)
density=sw_pden(salt,temp,abs(depth),0);

%momentum flux (Pa)
Fm=w.*w.*density;   

experiment.uTot=uTot;
experiment.melt=melt;
%% Put plume properties in output structure
experiment.plume.depth = depth;
experiment.plume.radius = radius;
experiment.plume.w = w; 
experiment.plume.temp = temp;
experiment.plume.salt = salt; 
experiment.plume.density = density;
experiment.plume.volumeFlux = Fv;
experiment.plume.momentumFlux = Fm; 
experiment.plume.melt = melt;
experiment.plume.type = options.type;

%find depth of maximum melt
[~,mi] = max(experiment.plume.melt);
experiment.plume.maximumMD = depth(mi); 

if strcmp(options.type,'line') || strcmp(options.type,'point')
    % find terminal level based on density
    tl = find(experiment.plume.density <= ambientRho,1);
    experiment.plume.neutralDensity = depth(tl);

    %find where integration stops
    wIndex = find(isnan(experiment.plume.w),1,'last');

    %if integration did not stop, then plume reached surface
    if isempty(wIndex)
        wIndex = 1;
    end

    %find maximum height (which is above terminal depth)
    experiment.plume.maximumH = depth(wIndex);

elseif strcmp(options.type,'ambient')
    % the above values have already been calculated separately for each plume
    experiment.plume.neutralDensity = nD;
    experiment.plume.maximumH = maxH;
end

%add units for all plume variables
experiment.plume.units={'depth (m)'; 'radius (m)'; 'w (m/s)'; 'temp (C)'; 'salt (psu)'; 'density (kg/m^3)'; 'volumeFlux (m^3/s)'; 'momentumFlux (Pa)'; 'melt (m/day)'; 'maximumMD (m)'; 'neutralDensity (m)'; 'maximumH (m)'};

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
    