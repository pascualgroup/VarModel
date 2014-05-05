% Get JSONlab for loadjson and savejson here:
% http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave
% addpath('path/to/jsonlab');

% Load existing parameters file into Matlab struct
params = loadjson('example_parameters.json');

% OR: create a new parameters struct from scratch
params = {};

% Modify parameters
params.dbFilename = 'output.sqlite'
params.withinHost = {};
params.withinHost.clearanceRateInitial = 0.5;

% Write struct to output JSON file
savejson('', params, 'filename-modified.json');







% {
% 	"dbFilename": "output.sqlite",
% 	"randomSeed": 0,
% 	"tYear": 365,
% 	"tEnd": 365,
% 	"seasonalUpdateEvery": 1,
% 	"genePoolSize": 30,
% 	"genesPerStrain": 5,
% 	"pRecombinant": 0,
% 	"tLiverStage": 7,
% 	"distanceFunction": {
% 		"power": 1
% 	},
% 	"hostLifetimeDistribution": {
% 		"dt": 365,
% 		"pdf": [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
% 	},
% 	"transmission": {
% 		"coinfectionReducesTransmission": 1
% 	},
% 	"withinHost": {
% 		"clearanceRateInitial": 0.5,
% 		"clearanceRateMidCourse": 0.01,
% 		"activationRate": 100,
% 		"deactivationRateImmune": 100,
% 		"deactivationRateNotImmune": 0.02
% 	},
% 	"populations": {
% 		"size": 1000,
% 		"nInitialInfections": 5,
% 		"bitingRate": {
% 			"mean": 1,
% 			"amplitude": 0
% 		},
% 		"immigrationRate": 0.01,
% 		"x": 0,
% 		"y": 0,
% 		"selfDistance": 0.5
% 	},
% 	"genes": {
% 		"transmissibility": 0.8,
% 		"immunityLossRate": 0.001,
% 		"clinicalImmunityLossRate": 0.1
% 	},
% 	"trackClinicalImmunity": 1
% }

