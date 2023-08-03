# BPTModel

## Description
BPTModel is a one-dimensional model based on buoyant plume theory in the context of marine-terminating glaciers, solved in depth space. This versatile model includes configurations for point-, line-, and stacked plumes featuring an option to include horizontal along-ice velocities and easy customization of a variety of coefficients and initial conditions.

The most current version of BPTModel can be found at: 

https://github.com/PolarCoasts/BPTModel

## Requirements
- MATLAB version 2019b or later
- seawater toolbox


## Quick start guide
The model can be run for a subglacial discharge plume with just a profile of ambient ocean properties (aDepth,aTemp,aSalt), a grounding line depth, and an estimate of subglacial discharge flux. All depth values should be entered as negative values, increasing towards the free surface. If the ambient profile includes data from near the free surface down to at least the grounding line depth, it can be input as is and the model, by default, will extrapolate to the surface and interpolate to a smaller depth interval. If there are large gaps in the data or the profile doesn't extend to the grounding line, the model defaults may lead to unexpected results (see "Inputs, outputs, and customization" for more information). The default configuration is a line plume with a 100 m wide discharge outlet. 

	exp1 = BPTmodel(aDepth,aTemp,aSalt,-100,300);
	
Model output (exp1) will be a structure with four substructures. The plume results can be found under exp1.plume, while the other substructures return the associated ambient ocean profile, coefficient values used, and initial conditions. The BPTModel package includes a plotting function for quick visualization of the model output. 

To plot plume radius, vertical velocity, and melt rate:

	[fig,ax] = BPTplots(exp1);
	
To plot all eight variables in the model results:

	[fig,ax] = BPTplots(exp1,'all');

Many features of the model are simple to customize.

To increase discharge outlet width to 200 m:
	
	exp2 = BPTmodel(aDepth,aTemp,aSalt,-100,300,W=200);
	
To run with point-plume geometry:

	exp3 = BPTmodel(aDepth,aTemp,aSalt,-100,300,type='point');
	
The included script, BPTexample, contains a sample ambient profile and demonstrates more of the model's features. Continue reading for a more thorough explanation of BPTModel, its capabilities and limitations. 


## Contents of package
The BPTModel package contains the following functions and scripts:
- **BPTmodel**:  model framework; accepts user inputs and calls underlying functions
- **line_plume**:  called by BPTmodel; solves system of ODEs for line-plume geometry (not user-callable)
- **point_plume**:  called by BPTmodel; solves system of ODEs for point-plume geometry (not user-callable)
- **melt_calc**:  called by BPTmodel, line_plume, and point_plume; solves for submarine melt rate associated with upwelling plume and can also be used as a stand-alone function to estimate submarine melt, neglecting plume dynamics
- **BPTplots**:	 plotting function; quick visualization of model output
- **BPTexample**:  contains sample ambient profile and demonstrates usage of several model features


## How the model works
### The theory...
Buoyant plume theory describes the evolution of plume properties with height. A melt parameterization describes the boundary layer physics and relates the plume properties to the melt rate. BPTmodel couples buoyant plume theory with a melt parameterization so that the plume properties control the melt rate and the melt feeds back into the plume as a source of buoyancy. Thus, the results represent the evolution of plume properties with height based on the entrainment of ambient water and the addition of melt water from the ice face. As such, the model is sensitive to the prescription of plume geometry.  Line plumes are modeled as wedge-shaped plumes emanating from a wide discharge outlet. Entrainment of ambient water is accounted for along the ice-parallel side of the plume, neglecting entrainment at the edges. Point plumes are modeled as half-conical plumes emanating from a narrow channel.

BPTmodel utilizes equations for the conservation of mass, momentum, heat, and salt, as described in Jenkins (2011) for line-plume geometry and modifications for point-plume geometry as described in Cowton et al. (2015). BPTmodel does not account for the slope of the ice face, so it is applicable to glaciers with near-vertical ice faces or undercut glaciers with a) slopes steep enough that the effects of friction are small and b) subglacial discharge fluxes high enough that the buoyancy flux from ice melt is small. In the case of undercut glaciers, the results represent the solution in depth coordinates, not along-ice coordinates.

### In more practical terms...
BPTmodel takes the input ambient profiles*, grounding line depth, and subglacial discharge flux and calculates initial plume radius and velocity according to default or user-defined settings. Initial temperature is set to the freshwater freezing point at the grounding line depth and initial salinity is set to near zero (0.0001 psu). These initial conditions as well as the boundary conditions (ambient profile and default or user-defined ice temperature) are then used to solve the conservation equations coupled with the melt parameterization by integrating from the grounding line to the surface (or to the point where momentum is zero). 

The resulting plume profile is then interpolated back onto a uniform depth interval (the ODE solver uses a variable depth interval) and additional variables are calculated. The ambient ocean profile used for the integration, the initial conditions, and all coefficients (default or user-defined), as well as the plume profile results are loaded as substructures in the output structure.

*By default, ambient profiles are linearly interpolated to 0.1 m depth intervals from the free surface to the grounding line and missing data at top and bottom of profile will be linearly extrapolated. Some of this behavior can be modified (see "Inputs, outputs, and customization").


## How to use BPTmodel
### Subglacial discharge plumes
In the simplest case, a user begins with a profile of ambient ocean conditions for the full water column, the glacier grounding line depth at the location of the discharge outlet, a value for subglacial discharge, and a preference for plume geometry. The default model behavior is designed to accommodate ambient ocean profiles that don't extend quite to the surface and may contain small gaps in data by linearly interpolating and extrapolating to the surface and grounding line. Some modifications can be made to this behavior (see "Inputs, outputs, and customization"), but any significant deficiencies in the ambient ocean profile should be handled by the user prior to input to the model. 

An ambient ocean profile (aDepth,aTemp,aSalt), grounding line depth, and subglacial discharge flux are the minimum inputs required. All depths should be entered as negative values, increasing towards the free surface. By default, the model will run as a line plume with a discharge outlet width of 100 m.

	exp1 = BPTmodel(aDepth,aTemp,aSalt,-100,300);
	
The outlet width (W) can be set with one additional input:

	exp2 = BPTmodel(aDepth,aTemp,aSalt,-100,300,W=200);
	
Alternatively, the plume geometry (type) can be modified:

	exp3 = BPTmodel(aDepth,aTemp,aSalt,-100,300,type='point');
	
Any number of settings can be customized (see "Inputs, outputs, and customization") in any order in the input fields. For example, one could run a point plume with custom entrainment and drag coefficients:

	exp4 = BPTmodel(aDepth,aTemp,aSalt,-100,300,type='point',alpha=0.08,Cd=2e-3);
	
### Stacked meltwater plumes, in the absence of subglacial discharge
BPTmodel simulates submarine melting in the absence of subglacial discharge as a set of stacked meltwater plumes. In this case, submarine meltwater is the only source of buoyancy. The model begins with a meltwater plume at the grounding line. For each plume that loses momentum prior to reaching the surface, a new plume is started directly above it. A non-zero value for subglacial discharge flux must be entered in order to initiate a meltwater plume, and a value of 1e-10 is likely sufficient (see Jackson et al., 2020). Stacked plumes are calculated using the line plume module. Default or user-defined outlet width will be overridden and replaced with a value of 1 m so that results will be per meter of terminus width. 

	exp5 = BPTmodel(aDepth,aTemp,aSalt,-100,1e-10,type='stacked');
	
To account for along-ice velocities driven by the ambient fjord circulation, a horizontal velocity can also be included. The horizontal velocity will not impact the plume dynamics directly, but will be combined with the upwelling velocity to calculate the total along-ice velocity that drives melting (see Jackson et al., 2020).

	exp6 = BPTmodel(aDepth,aTemp,aSalt,-100,1e-10,type='stacked',uh=0.1);
	
### Plotting model results
A simple plotting function is included for quick visualization of the model output. When the only input is the results structure from the model, plume radius, vertical velocity, and melt rate will be plotted by default. The term 'all' can be used as a second input to get a figure with all available plume variables.

	[fig,ax] = BPTplots(exp1,'all');
	
Alternatively, the user may select specific variables:

	[fig,ax] = BPTplots(exp1,["w", "temp", "salt", "melt"]);
	
The variables available for plotting are:

| name			| description									| units	|
| :---				| :---											| :---:		|
| "radius"			| plume radius (point plume) or thickness (line plume) 	| m		|
| "w"				| plume vertical velocity							| m/s		|
| "temp"			| plume temperature*								| C		|
| "salt"			| plume salinity*									| psu		|
| "density"			| plume density*								| kg/m^3	|
| "volumeFlux"		| vertical flux of volume within plume					| m^3/s	|
| "momentumFlux"	| vertical flux of momentum within plume				| kg*m/s^2 |
| "melt"			| melt rate										| m/day	|
	
*These variables will be plotted with the ambient ocean profile. To omit the ambient ocean profile, enter 0 as the third input.

### Submarine melt rates independent of plume dynamics
The function that handles the melt parameterization for the model can also be used as a stand-alone function to estimate submarine melt rates neglecting plume dynamics. This function can handle a profile of horizontal velocities and includes several customizable variables (see "Inputs, outputs, and customization").

	melt=melt_calc(aDepth,uprofile,aTemp,aSalt,iceT=0);


## Inputs, outputs, and customization
### Inputs to BPTmodel 	
| name (*required)	| description						| details 										| units 		| default	|
| :---				| :---								| :---											| :---:			| :---:		|
| aD*				| ambient ocean depth profile			| negative values, increasing towards surface 			| m			|		|
| aT*				| ambient ocean temperature profile		|											| C			|		|
| aS*				| ambient ocean salinity profile			|											| psu			|		|
| GL*				| grounding line depth				| negative value 								| m			|		|
| Q*				| subglacial discharge flux				|											| m^3/s		|		|
| type			| plume geometry/configuration			|											|			| 'line'	|
| iceT			| ice temperature					|											| C			| -10		|
| W				| discharge outlet width				| for line plumes 								| m			| 100		|
| uh				| ambient ocean horizontal velocity		|											| m/s			| 0		|
| ui				| initial plume velocity				|											| m/s			| 'balance'	|
| alpha			| entrainment coefficient				|											|			| 0.1		|
| Cd				| drag coefficient					|											|			| 2.5e-3	|
| gammaT			| thermal transfer coefficient			|											|			| 0.022	|
| gammaS			| haline transfer coefficient				|											|			| 0.00062	|
| intMethod		| interpolation method				|											|			| 'linear'	|
| extValue_T		| temperature extrapolation value		|											|			| 'extrap'	|
| extValue_S		| salinity extrapolation value			|											|			| 'extrap'	|

### Customizing inputs to BPTmodel
Any optional inputs can be modified from the default by entering 'name=value' as an input argument into the BPTmodel. Below is a brief description of each and its role in the model framework. See table above for default values.
	
- type:  plume geometry/configuration ('line', 'point', or 'stacked')
    - BPTmodel can be run with either a 'line' or 'point' plume geometry.  The line plume geometry assumes that subglacial discharge is distributed uniformly across a specified portion of the grounding line (discharge outlet width). This configuration only accounts for entrainment along the width (ice-parallel dimension) of the plume, neglecting entrainment at the ends. In this case, the 'radius' of the plume refers to its thickness (ice-perpendicular dimension). The point plume geometry assumes that subglacial discharge comes from a localized channel.
    - Setting type='stacked' will run the model as a series of stacked line plumes. Use this option to simulate meltwater plumes in the absence of subglacial discharge. Note that a non-zero value for subglacial discharge is required to initiate a plume (ex. 1e-10). Discharge outlet width will be set to 1 m so that results will be per meter of terminus width.

- iceT:  ice temperature			
    - In the melt parameterization, heat is consumed to bring ice up to the freezing point before the phase change occurs. The model will accept any value less than or equal to zero.
	
- W:  discharge outlet width	
    - In a line plume, the discharge outlet width is the along-terminus distance over which the subglacial discharge is uniformly distributed. The model calculates plume properties per unit width. Model results are then multiplied by the outlet width for the volumetric properties of volume and momentum fluxes. Stacked melt plumes are always calculated per m terminus width, so the default or user input for outlet width is overwritten with a value of 1. Outlet width is not applicable for point plumes, so the default or user input for outlet width is ignored.

- uh:  ambient ocean horizontal velocity
    - This is the horizontal along-ice velocity that will be combined with the vertical plume velocity to get a total along-ice velocity for the melt parameterization. The horizontal velocity is not accounted for in the plume model, only the coupled melt parameterization. There is no direct affect on the plume dynamics, but increased melt rates due to the inclusion of a horizontal velocity will increase the buoyancy flux to the plume. 

- ui:  initial plume velocity (scalar value or 'balance')
    - This is the initial vertical velocity of the plume. The user may input a numerical value or use the default option 'balance', where the initial velocity is calculated as a balance of momentum and buoyancy:
        - For a line plume: $u_i= {{g'q} \over \alpha}$
        - For a point plume*: $u_i= {{2 \over \pi} ({{ \pi^2 g'} \over 8\alpha})^{2 \over 5} Q^{1 \over 5}}$
    - The model assumes all velocity is upward so that the initial radius of the plume can be calculated by mass conservation for the subglacial discharge flux at the initial velocity.
*The initial velocity for a point plume based on a balance of momentum and buoyancy comes directly from the conservation equations as applied in this model and is slightly different from that which appears in the literature (e.g., Slater et. al, 2016). This results in a reduction in initial velocity of less than 10\%. Model results away from the grounding line are less sensitive to this choice.

- alpha:  entrainment coefficient
    - This controls the rate at which ambient ocean water is entrained into the plume.
	
- Cd:  drag coefficient
    - Drag coefficient influences both the loss of momentum to the ice wall and the exchange of heat and salt across the ice-ocean boundary.
	
- gammaT:  thermal transfer coefficient
    - This controls the rate of heat transfer across the ice-ocean boundary
	
- gammaS:  haline transfer coefficient
    - This controls the rate of salt transfer across the ice-ocean boundary
	
- intMethod:  interpolation method*
    - The user input ambient ocean profiles will be interpolated to 0.1 m depth intervals from the free surface to the grounding line using the built-in MATLAB function interp1. Any interpolation method available to interp1 can be used here. 
	
- extValue_T:  temperature extrapolation value* (scalar value or 'extrap')
    - If the range of the input ambient ocean profile does not extend to the free surface and/or the grounding line, the data will be extrapolated as per interp1 unless a constant value is given. If a constant value is given, the same value will be used for missing data at the top and bottom of the water column.
	
- extValue_S:  salinity extrapolation value* (scalar value or 'extrap')
    - If the range of the input ambient ocean profile does not extend to the free surface and/or the grounding line, the data will be extrapolated as per interp1 unless a constant value is given. If a constant value is given, the same value will be used for missing data at the top and bottom of the water column.
	
*The default methods for interpolation and extrapolation of the ambient ocean profile work well for profiles with good coverage from near the surface down to, and including, the grounding line depth. For profiles with well-understood deficiencies, the user may opt to handle them within the model by modifying the interpolation and extrapolation methods. To avoid unexpected behavior, it is recommended that ambient ocean profiles that do not extend to the grounding line or that have significant gaps in data be handled prior to being input to the model.

### Outputs from BPTmodel
The model output is packaged into a single structure with the following four substructures.

<br/>

**plume**:  Substructure containing the model results.
| field name		| description												| units		| calculation		|
| :---				| :---														| :---:			| :---:				|
| depth			| depth coordinate of all profiles									| m			|				|
| radius			| plume radius (for line plumes, this is the ice-perpendicular thickness)	| m			|				|
| w				| vertical velocity											| m/s			|				|
| temp			| plume temperature											| C			|				|
| salt				| plume salinity												| psu			|				|
| density			| plume density												| kg/m^3		|				|
| area			| cross-sectional area of plume									| m^2		|				|
| volumeFlux		| volume of flow through cross-sectional area of plume per unit time		| m^3/s		| $F_V=Aw$		|
| momentumFlux	| flow of momentum through cross-sectional area of plume per unit time	| kg m/s^2		| $F_m=ρAw^2$	|
| melt			| ice melt rate												| m/day		|				|
| plumeType		| geometry/configuration ('line', 'point', or 'stacked')					|			|				|
| maximumMD		| depth of maximum melt rate									| m			|				|
| neutralDensity		| depth at which plume density was equal to that of ambient ocean		| m			|				|
| maximumH		| shallowest depth (maximum height) reached by plume				| m			|				|
| units			| units for each of the plume variables above						|			|				|

<br/>

**ambient_ocean**:  Substructure containing the ambient ocean profiles used by the model. This is the interpolated/extrapolated version of the user-input profile variables and additional variables calculated from them.
| field name		| description												| units		|
| :---				| :---														| :---:			|
| depth			| depth coordinate of all profiles									| m			|
| temp			| ambient ocean temperature									| C			|
| salt				| ambient ocean salinity										| psu			|
| density			| ambient ocean density										| kg/m^3		|
| u_horiz			| horizontal along-ice velocity									| m/s			|
| N2				| buoyancy frequency										| s^-1		|
| units			| units for each of the ambient ocean variables above				|			|

<br/>

**initial_cond**:  Substructure containing the user-input or default variables for initial conditions used by the model.
| field name		| description												| units		|
| :---				| :---														| :---:			|
| GL				| grounding line depth 										| m			|
| outletW			| discharge outlet width (line plume)								| m			|
| Q				| subglacial discharge flux										| m^3/s		|
| iceTemp			| temperature of ice											| C			|
| ti				| initial plume/discharge temperature								| C			|
| si				| initial plume/discharge salinity									| psu			|
| ui				| initial plume vertical velocity or method of calculation*				| m/s 		|
| plumeType		| geometry/configuration ('line', 'point', or 'stacked')					|			|
| units			| units for each of the initial condition variables above				|			|

*If this field reads 'balance', then the initial velocity was calculated as a balance of momentum and buoyancy

<br/>

**coeff**:  Substructure containing the user-input or default values for coefficients used by the model.
| field name		| description				|
| :---				| :---						|
| alpha			| entrainment coefficient		|
| Cd				| drag coefficient			|
| gammaT			| thermal transfer coefficient	|
| gammaS			| haline transfer coefficient		|

### Inputs to melt_calc
| name (*required)	| description						| details								| units	| default	|
| :---				| :---								| :---									| :---:		| :---:		|
| z*				| depth coordinate					| negative values, increasing towards surface	| m		|		|
| u*				| along-ice velocity that drives melting	| combination of vertical and horizontal velocity	| m/s		|		|
| T*				| ambient temperature				| 									| C		|		|
| S*				| ambient salinity					|									| psu		|		|
| iceT			| ice temperature					|									| C		| -10		|
| Cd				| drag coefficient					|									|		| 2.5e-3	|
| gammaT			| thermal transfer coefficient			|									|		| 0.022	|
| gammaS			| haline transfer coefficient				|									|		| 0.00062	|
| cw				| heat capacity of seawater			|									| J/kg/C	| 3974	|
| ci				| heat capacity of ice					|									| J/kg/C	| 2009	|
| L				| latent heat of fusion					|									| J/kg	| 335000	|
| lambda1			| variation of freezing point with salinity	|									| C/psu	| -0.0573	|
| lambda2			| freezing point offset					|									| C		| 0.0832	|
| lambda3			| variation of freezing point with depth	|									| C/m	| 7.61e-4	|

### Outputs from melt_calc

	melt=melt_calc(z,u,T,S);
	
	[melt,Tb,Sb]=melt_calc(z,u,T,S);
	
**melt**:  Ice melt rate in m/s

**Tb**: temperature at ice-ocean boundary (salinity- and depth-dependent freezing point)

**Sb**: salinity at ice-ocean boundary (dependent upon Tb and melt)


## Future developments
- Update calculations of thermodynamic properties with GSW toolbox (McDougall and Barker, 2011)
- Modify BPTmodel to allow for profile of horizontal velocity
- Account for terminus slope in conservation equations
- Add modules for alternate plume geometries
- Develop a version of the model for Python


## Credits
The BPTModel package was developed by Bridget Ovall. The mathematical basis for this model comes from Jenkins (2011) with modifications for point-plume geometry from Cowton et. al (2015). Adaptations for stacked melt plumes with consideration of fjord-scale circulation are based on Jackson et. al (2020). 

If using this model in a published work, please provide a link to our lab's GitHub page so that others can find it, too. 

https://github.com/PolarCoasts/BPTModel

## References

Cowton, T., Slater, D., Sole, A., Goldberg, D., and Nienow, P. (2015). Modeling the impact of glacial runoff on fjord circulation and submarine melt rate using a new subgrid-scale parameterization for glacial plumes. Journal of Geophysical Research: Oceans, 120(2), 796–812. https://doi.org/10.1002/2014JC010324

Jackson, R. H., Nash, J. D., Kienholz, C., Sutherland, D. A., Amundson, J. M., Motyka, R. J., Winters, D., Skyllingstad, E., & Pettit, E. C. (2020). Meltwater Intrusions Reveal Mechanisms for Rapid Submarine Melt at a Tidewater Glacier. Geophysical Research Letters, 47(2). https://doi.org/10.1029/2019GL085335

Jenkins, A. (2011). Convection-Driven Melting near the Grounding Lines of Ice Shelves and Tidewater Glaciers. Journal of Physical Oceanography, 41(12), 2279–2294. https://doi.org/10.1175/JPO-D-11-03.1

McDougall, T.J., and P.M. Barker, 2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW) Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127, ISBN 978-0-646-55621-5.
