function experiment = BPTmodel(aD,aT,aS,GL,Q,options)

arguments
   aD (1,:) {mustBeLessThanOrEqual(aD,0)}                               % ambient profile depths (m)
   aT (1,:) {mustBeEqualSize(aD,aT)}                                    % ambient profile of temperature (degrees C)
   aS (1,:) {mustBeEqualSize(aD,aS)}                                    % ambient profile of salinity (psu)
   GL (1,1) {mustBeLessThan(GL,0)}                                      % grounding line depth (m)
   Q (1,1) {mustBeGreaterThan(Q,0)}                                     % subglacial discharge flux (m^3/s)
   options.iceT (1,1) {mustBeLessThanOrEqual(options.iceT,0)}=-10       % ice temperature (degrees C)
   options.W (1,1) {mustBeGreaterThanOrEqual(options.W,1)}=100          % line plume outlet width (m) 
   options.uh (1,1) double=0                                            % horizontal velocity relevant for melting (m/s)
   options.ui {mustBeValueOrChar(options.ui,"balance")}='balance'       % initial velocity (m/s), 'balance' of momentum and buoyancy or set value
   options.alpha (1,1) double=0.1                                       % entrainment coefficient
   options.Cd (1,1) double=2.5e-3                                       % drag coefficient
   options.gammaT (1,1) double=0.022                                    % thermal transfer coefficient
   options.gammaS (1,1) double=0.00062                                  % haline transfer coefficient
   options.intMethod char {mustBeMember(options.intMethod,{'linear','nearest','next','previous','pchip','cubic','v5cubic','makima','spline'})}='linear'  % method for interpolating ambient profile
   options.extValue_T {mustBeValueOrChar(options.extValue_T,"extrap")}='extrap' % to extend ambient profile to full plume range, extrapolate or set to specific value
   options.extValue_S {mustBeValueOrChar(options.extValue_S,"extrap")}='extrap' % to extend ambient profile to full plume range, extrapolate or set to specific value
   options.type char {mustBeMember(options.type,{'line','point','stacked','corner'})}='line' % run model as line- or point-plume, stacked runs a series of line-plumes to simulate melt with no discharge
    %NOTE: undocumented type='corner' is experimental
end
% BPTMODEL 1D model of buoyant plumes in the context of marine-terminating glaciers, solved in depth space
%
%   experiment = BPTmodel(aD,aT,aS,GL,Q) for row vectors aD, aT, and aS (profiles of ambient ocean depth, temperature, and salinity)
%   and scalars GL (grounding line) and Q (subglacial discharge flux). All depth values (aD and GL) should be given as negative
%   increasing towards the free surface. Q should be greater than zero (use a small value (eg. 1e-10) for 'stacked' melt plumes). 
%   Output structure contains four substructures that contain model results and all inputs.
% 
%   experiment = BPTmodel(aD,aT,aS,GL,Q,opt) Specify values for optional inputs shown below (default):     
%     iceT (-10)                  temperature of ice (C)
%     W (100)                     discharge outlet width for line plumes (m)
%     uh (0)                      horizontal velocity to be used for calculating melt (m/s)
%     ui ('balance')              initial upwelling velocity (specify value or 'balance')
%     alpha (0.1)                 entrainment coefficient
%     Cd (2.5e-3)                 drag coefficient
%     gammaT (0.022)              thermal transfer coefficient
%     gammaS (0.00062)            haline transfer coefficient
%     intMethod ('linear')        interpolation method for ambient profile (methods of interpolation for interp1)
%     extValue_T ('extrap')       value to fill ambient temperature profile with to extend to full water depth (or specify 'extrap')
%     extValue_S ('extrap')       value to fill ambient salinity profile with to extend to full water depth (or specify 'extrap')
% 
% Examples: 
%   % minimum inputs
%   experiment = BPTmodel(aD,aT,aS,GL,Q)  
%
%   % setting two optional values, run as point plume
%   experiment = BPTmodel(aD,aT,aS,GL,Q,type='point',Cd=2e-3,alpha=0.11) 
% 
% OUTPUT (structure with four substructures):
%     plume:
%     - depth/temp/salt/density     properties of upwelling plume
%     - w                           vertical velocity of plume
%     - melt                        melt rate
%     - neutralDensity              depth where plume density equals ambient density
%     - maximumH                    depth of maximum height reached by plume (momentum=0 or surface)
%     - radius                      thickness of plume (perpendicular to ice)
%     - area                        area that plume occupies at each depth (radius*outlet width or pi/2*radius^2)
%     - maximumMD                   depth of maximum melt rate
%     - volumeFlux                  volume flux of plume
%     - momentumFlux                momentum flux of plume
%     - plumeType                   plume geometry used (ex., 'point' or 'line')
%     - units                       list of units for plume variables
% 
%     ambient_ocean:
%     - depth/temp/salt             input profiles interpolated to 0.1 m increments and extended to full water column
%     - u_horiz                     input horizontal velocity
%     - density                     density of ambient profile
%     - N2                          buoyancy frequency of ambient profile
%     - units                       list of units for ambient variables
%        
%     initial_cond:
%     - GL                          grounding line depth
%     - outletW                     discharge outlet width
%     - Q                           subglacial discharge flux
%     - iceTemp                     temperature of ice
%     - ti                          initial plume temperature
%     - si                          initial plume salinity
%     - ui                          initial velocity at grounding line depth
%     - plumeType                   plume geometry used (ex., 'point' or 'line')
%     - units                       list of units for glacier variables
%     
%     coeff:
%     - alpha                       entrainment coefficient
%     - Cd                          drag coefficient
%     - gammaT                      thermal transfer coefficient
%     - gammaS                      haline transfer coefficient
% 
% FUNCTIONS called:
%       line_plume(...) or point_plume(...)
%       melt_calc(...)
%
%
% For the most current version, see <a href="matlab: 
% web('https://github.com/BridgetOvall/BPTModel.git')">BPTmodel on GitHub</a>.
%

% supress integration failure warning that occurs for any plume that does not reach the surface
warning('off','MATLAB:ode45:IntegrationTolNotMet')

% for stacked plumes, always use a unit width
if strcmp(options.type,'stacked')
    if options.W~=100 && options.W~=1
        warning('Stacked plumes are always calculated per meter terminus width')
    end
    options.W=1;
end

% convert discharge to m2/s for line plume or keep as m3/s for point plume
if strcmp(options.type,'line') || strcmp(options.type,'stacked')
    qsg = Q/options.W;
elseif strcmp(options.type,'point') || strcmp(options.type,'corner')
    qsg = Q;
end

% simple validations of input profiles
% --- This is designed to catch the most obvious problems, but should not be used as a validation of adequate data ---
top=max(aD(~isnan(aT) & ~isnan(aS))); bottom=min(aD(~isnan(aT) & ~isnan(aS))); getinput=0;
if sum(~isnan(aT))<10 || sum(~isnan(aS))<10
    error(sprintf("Ambient ocean profiles must contain more than 10 good data points.\n"))
end
if isstring(options.extValue_T) || ischar(options.extValue_T)
    TintMsg="Missing temperature data will be extrapolated per the methods set by interp1(). ";
else
    TintMsg="Missing temperature data will be filled with a value of "+options.extValue_T+". ";
end
if isstring(options.extValue_S) || ischar(options.extValue_S)
    SintMsg="Missing salinity data will be extrapolated per the methods set by interp1(). ";
else
    SintMsg="Missing salinity data will be filled with a value of "+options.extValue_S+". ";
end
if top<-10
    warning(sprintf("Ambient ocean profiles are missing data over the top "+ abs(top) + " m. " + TintMsg + SintMsg +" \n"))
    getinput=1;
end
if bottom-GL>5
    warning(sprintf("Ambient ocean profiles are missing data over the bottom " + (bottom-GL) + " m (above the grounding line). " + TintMsg + SintMsg +" \n"))
    getinput=1;
end
if sum(isnan(aT))/length(aT)>.1 || sum(isnan(aS))/length(aS)>.1
    warning(sprintf("Ambient ocean profiles contain greater than 10%% NaNs. \n"))
    getinput=1;
end
if abs(mean(diff(aD)))>5
    warning(sprintf("Ambient ocean profiles have a mean depth interval greater than 5 m\n"))
    getinput=1;
end
if getinput==1
    userinfo=input("Would you like to continue? (y/n)",'s');
    if isempty(userinfo)
        userinfo='n';
    end
    if userinfo~='y'
        experiment=[];
        fprintf("\nModel run aborted\n")
        return
    end
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
experiment.ambient_ocean.depth = depth;
experiment.ambient_ocean.temp = ambientTemp;
experiment.ambient_ocean.salt = ambientSalt;
experiment.ambient_ocean.density = ambientRho;
experiment.ambient_ocean.u_horiz = options.uh;

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

% find ambient rho and plume rho at grounding line --> calculate g' (gp)
iGL = find(abs(depth)==abs(GL),1);
ambientRho_GL = ambientRho(iGL);
plumeRho_GL = sw_pden(si,ti,abs(GL),0);
gp = (ambientRho_GL-plumeRho_GL)*const.g/const.rho0;

% initial velocity of plume
if strcmp(options.type,'line') || strcmp(options.type,'stacked')
    if ~isscalar(options.ui) && ismember(options.ui,"balance")
        ui = (gp*qsg/options.alpha)^(1/3); 
    else
        ui = options.ui;
    end
    bi = qsg / ui; %plume thickness in cross-terminus direction (m)
elseif strcmp(options.type,'point')
    if ~isscalar(options.ui) && ismember(options.ui,"balance")
        ui = 2/pi*(pi^2*gp/(8*options.alpha))^(2/5)*qsg^(1/5); 
    else
        ui = options.ui;
    end
    bi = sqrt(2*qsg / (pi*ui)); %plume radius (m)
elseif strcmp(options.type,'corner')
    if ~isscalar(options.ui) && ismember(options.ui,"balance")
        ui = 4/pi*(pi^2*gp/(32*options.alpha))^(2/5)*qsg^(1/5);
    else
        ui = options.ui;
    end
    bi = sqrt(4*qsg / (pi*ui)); %plume radius (m)
end
%% Solve coupled system of ODEs 
% using a 4th order MATLAB integrator
%    z = depth
%    X is m x 4 with columns of plume: width, velocity, T , S

if strcmp(options.type,'line')
    [z,X]=ode45(@(Z,x) line_plume(Z,x,const,options,experiment.ambient_ocean),[GL,surface],[bi,ui,ti,si]);
elseif strcmp(options.type,'point')
    [z,X]=ode45(@(Z,x) point_plume(Z,x,const,options,experiment.ambient_ocean),[GL,surface],[bi,ui,ti,si]);
elseif strcmp(options.type,'corner')
    [z,X]=ode45(@(Z,x) corner_plume(Z,x,const,options,experiment.ambient_ocean),[GL,surface],[bi,ui,ti,si]);
elseif strcmp(options.type,'stacked')
    di=GL; % initial depth of bottom plume is grounding line depth
    count=0; maxH=[]; nD=[];
    while di<0
        count=count+1;
        [tz,tX]=ode45(@(Z,x) line_plume(Z,x,const,options,experiment.ambient_ocean),[di,surface],[bi,ui,ti,si]);
        [tz,ia]=unique(tz); %occassionally, intervals are too small to be distinguishable
        tX=tX(ia,:);
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
            fprintf('    Solving for next plume ... \n')
        else
            fprintf('    Plume integration complete. Preparing final output...\n')
        end
    end
end

%% Interpolate to specified grid 
if strcmp(options.type,'stacked')
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
    radius=sum(ar,1,'omitnan');
    w=sum(aw,1,'omitnan');
    temp=sum(at,1,'omitnan');
    salt=sum(as,1,'omitnan');
else
    [z,ia]=unique(z); %occassionally, intervals are too small to be distinguishable
    X=(X(ia,:));
    radius=interp1(z,X(:,1),depth);
    w=interp1(z,X(:,2),depth);
    temp=interp1(z,X(:,3),depth);
    salt=interp1(z,X(:,4),depth);
end

%% Compute additional variables from solution of ODEs

% define a geometric factor to account for differing plume areas
if strcmp(options.type,'line') || strcmp(options.type,'stacked')
    geom = options.W;
elseif strcmp(options.type,'point')
    geom = pi.*radius./2;
elseif strcmp(options.type,'corner')
    geom = pi.*radius./4;
end

% area (m^2)
area = geom.*radius;

% volume flux (m^3/s)
Fv = area.*w;

% calculate total velocity that drives melting: upwelling + horizontal
uTot=sqrt(w.^2 + options.uh.^2);

% melt rate (m/s --> m/day)
[melt_ms,~,~] = melt_calc(depth, uTot, temp, salt, lambda1=const.lambda1, lambda2=const.lambda2, lambda3=const.lambda3, gammaS=options.gammaS, gammaT=options.gammaT, Cd=options.Cd, iceT=options.iceT, L=const.L, ci=const.ci, cw=const.cw);
melt=24*60*60*melt_ms;

%plume density (kg m^-3)
density=sw_pden(salt,temp,abs(depth),0);

%momentum flux (kg m/s^2)
Fm=w.*w.*density.*area;   

%% Put plume properties in output structure
experiment.plume.depth = depth;
experiment.plume.radius = radius;
experiment.plume.w = w; 
experiment.plume.temp = temp;
experiment.plume.salt = salt; 
experiment.plume.density = density;
experiment.plume.area = area;
experiment.plume.volumeFlux = Fv;
experiment.plume.momentumFlux = Fm; 
experiment.plume.melt = melt;
experiment.plume.plumeType = options.type;

%find depth of maximum melt
[~,mi] = max(experiment.plume.melt);
experiment.plume.maximumMD = depth(mi); 

if strcmp(options.type,'line') || strcmp(options.type,'point') || strcmp(options.type,'corner')
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

elseif strcmp(options.type,'stacked')
    % the above values have already been calculated separately for each plume
    experiment.plume.neutralDensity = nD;
    experiment.plume.maximumH = maxH;
end

%add units for all plume variables
experiment.plume.units={'depth (m)'; 'radius (m)'; 'w (m s^{-1})'; 'temp (C)'; 'salt (psu)'; 'density (kg m^{-3})'; 'area (m^2)'; 'volumeFlux (m^3 s^{-1})'; 'momentumFlux (kg m s^{-2})'; 'melt (m day^{-1})'; 'maximumMD (m)'; 'neutralDensity (m)'; 'maximumH (m)'};

%% calculate N2 of ambient profile
experiment.ambient_ocean.N2 = real(sqrt(const.g/const.rho0 * diff(ambientRho))); 

%add units for all ambient ocean variables
experiment.ambient_ocean.units={'depth (m)'; 'temp (C)'; 'salt (psu)'; 'density (kg m^{-3})'; 'u_horiz (m s^{-1})'; 'N2 (s^{-1})'};

%% Return initial conditions
experiment.initial_cond.GL = GL;
if strcmp(options.type,'line') || strcmp(options.type,'stacked')
    experiment.initial_cond.outletW = options.W;
elseif strcmp(options.type,'point') || strcmp(options.type,'corner')
    experiment.initial_cond.outletW = 'not applicable';
end    
experiment.initial_cond.Q = Q;
experiment.initial_cond.iceTemp = options.iceT;
experiment.initial_cond.ti = ti;
experiment.initial_cond.si = si;
experiment.initial_cond.ui = options.ui;
experiment.initial_cond.plumeType = options.type;
experiment.initial_cond.units={'GL (m)'; 'outletW (m)'; 'Q (m^3 s^{-1})'; 'iceTemp (C)'; 'ti (C)'; 'si (psu)'; 'ui (m s^{-1})'};

%% Return coefficients used
experiment.coeff.alpha = options.alpha;
experiment.coeff.Cd = options.Cd;
experiment.coeff.gammaT = options.gammaT;
experiment.coeff.gammaS = options.gammaS;

%%
% reinstate integration failure warning
warning('on','MATLAB:ode45:IntegrationTolNotMet')

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
    if ~isscalar(a)
        if ~ismember(a,string)
            eid='Value:Error';
            msg="extValue must be 'extrap' or a constant value";
            throwAsCaller(MException(eid,msg))
        end
    end
end
    