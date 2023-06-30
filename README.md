# BPTModel
 One dimensional model based on buoyant plume theory, applied to marine-terminating glaciers

## BPTmodel ( ) is a command-line function to run the one-dimensional model with a choice of line- or point-plume geometry
- it requires an input of an ambient profile (temp, salinity, depth), a grounding line depth, and a subglacial discharge flux
- alternatively, model can be run for ambient melt in the absence of subglacial discharge
	- 'ambient' models melt as stacked line plumes, with each plume originating directly above the detrainment depth of the previous
	- requires a minimal discharge flux to initiate the plume (ex: 1e-10 with width=1m)
- a variety of constants/coefficients are loaded with default values, but can be customized as optional inputs
- model results are output into a nested structure under sub-structure 'plume'
- values used to run the model (based on input or defaults) are also included in output under substructures 'ambient', 'glacier', and 'coeff'
	
### Examples:

	exp1=BPTmodel(aDepth,aTemp,aSalt,100,300);
		Model is run with ambient profile (aDepth,aTemp,aSalt) starting from the grounding line depth of 100 m with subglacial discharge flux of 300 m^3/s and defaults for all constants and coefficients
		
	exp2=BPTmodel(aDepth,aTemp,aSalt,100,300,W=150,alpha=0.08)
		Model is run as in previous example, but with outlet width set at 150 m and entrainment coefficient set at 0.08. All other constants and coefficients use default values
		
	exp3=BPTmodel(aDepth,aTemp,aSalt,100,300,type='point')
		Model is run as in exp1, but with point-plume geometry. Note that outlet width is not used in this geometry, so there is no need to override the default.
		
	exp4=BPTmodel(aDepth,aTemp,aSalt,100,1e-10,W=1,type='ambient')
		Model is run as a series of stacked line plumes to replicate ambient melt. Width is set to 1 to produce results per meter of terminus length and a very small discharge flux is input to initiate plumes.
		
## melt_calc(), used within BPTmodel(), can be used as a stand-alone function to estimate melt rates
- it requires an input of near-ice conditions (profiles of depth, along-ice velocity, ambient temperature and salinity) 
	- will also accept scalar values of each
- a variety of constants/coefficients are loaded with default values, but can be customized as optional inputs
- outputs melt rate in m/s (optional outputs of boundary T and S)
		
		