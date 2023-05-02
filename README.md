# BPTModel
 One dimensional model based on buoyant plume theory, applied to marine-terminating glaciers

## BPTmodel ( ) is a command-line function to run the one-dimensional model with a line-plume geometry
- it requires an input of an ambient profile (temp, salinity, depth), a grounding line depth, and a subglacial discharge flux
- a variety of constants/coefficients are loaded with default values, but can be customized as optional inputs
- model results are output into a nested structure under sub-structure 'plume'
- values used to run the model (based on input or defaults) are also included in output under substructures 'ambient', 'glacier', and 'coeff'
	
### Examples:

	exp1=BPTmodel(aDepth,aTemp,aSalt,100,300);
		Model is run with ambient profile (aDepth,aTemp,aSalt) starting from the grounding line depth of 100 m with subglacial discharge flux of 300 m^3/s and defaults for all constants and coefficients
		
	exp2=BPTmodel(aDepth,aTemp,aSalt,100,300,W=150,alpha=0.08)
		Model is run as in previous example, but with outlet width set at 150 m and entrainment coefficient set at 0.08. All other constants and coefficients use default values
		
		