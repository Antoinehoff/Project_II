#!/bin/bash



if [ ! -e "output" ]; then
  mkdir output
fi
#Parameters:
appearanceNormWeight=100.0 #not very influent
complianceMaxFactor=1.2
exemplarDownsampling=2 #Reduce the size of the pattern input w.r.t. the domain
exponentSimilarityMetric=2.6 #higher -> force the ressemblance to the pattern
filterRadius=1.1 #smaller -> higher grain definition of the output
neighborhoodSize=20 #smaller -> will try to replicate more the pattern given
patchMatchIter=50 #higher -> more pattern accurate
penalty=3.0 #
treshPedersen=0.3 # No clue of what it does

name="bridge_non_symmetric_HD"
img_name="${name}.png"
txt_name="${name}.txt"
json_name="${name}.json"
echo ${outputname}
  echo "{
	\"appearanceNormWeight\": ${appearanceNormWeight},
	\"complianceMaxFactor\": ${complianceMaxFactor},
	\"densityMin\": 0.0,
	\"exemplarDownsampling\": ${exemplarDownsampling},
	\"exemplarPath\": \"../img/grid.ppm\",
	\"exponentSimilarityMetric\": ${exponentSimilarityMetric},
	\"filterRadius\": ${filterRadius},
	\"filterType\": \"FilterType.Density\",
	\"hasGui\": false,
	\"interpolationType\": \"InterpolationType.SIMP\",
	\"lambdaOccurrenceMap\": 20,
	\"maxSolverStep\": 100,
	\"neighborhoodSize\": ${neighborhoodSize},
	\"nelx\": 100,
	\"nely\": 50,
	\"numLevels\": 3,
	\"patchMatchIter\": ${patchMatchIter},
	\"penalty\": ${penalty},
	\"poissonCoeff\": 0.4,
	\"problemModule\": \"bridge\",
	\"problemType\": \"ProblemType.AppearanceWithMaxComplianceAndSymmetry\",
	\"treshPedersen\": ${treshPedersen},
	\"volumeFracMax\": 0.45,
	\"volumeFracMin\": 0,
	\"youngModulusMax\": 1.0,
	\"youngModulusMin\": 1e-09,
	\"materialDensity\": 0.0,
	\"a\":[[0,1]],
	\"c\":[[0.5,0.5]]
}" > input/${json_name}
python main.py input/${json_name} output/${img_name} > output/${txt_name}
