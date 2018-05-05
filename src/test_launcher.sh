#!/bin/bash



if [ ! -e "tests_folder" ]; then
  mkdir tests_folder
fi
#Parameters to study :
appearanceNormWeight=100.0 #not very influent
complianceMaxFactor=1.2
exemplarDownsampling=3.0 #Reduce the size of the pattern input w.r.t. the domain
exponentSimilarityMetric=1.2 #higher -> force the ressemblance to the pattern
filterRadius=1.1 #Control the smoother of the gradient in a local neighborhood
neighborhoodSize=20 #smaller -> will try to replicate more the pattern given
patchMatchIter=50 #higher -> more pattern accurate
penalty=4.0 #
treshPedersen=0.3 # No clue of what it does

array=(0.1 0.3 0.6 0.9)
for i in ${array[@]}; do
treshPedersen=${i}
outputname="treshPedersen${i}.png"
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
	\"nelx\": 50,
	\"nely\": 25,
	\"numLevels\": 3,
	\"patchMatchIter\": ${patchMatchIter},
	\"penalty\": ${penalty},
	\"poissonCoeff\": 0.4,
	\"problemModule\": \"bridge\",
	\"problemType\": \"ProblemType.AppearanceWithMaxCompliance\",
	\"treshPedersen\": ${treshPedersen},
	\"volumeFracMax\": 0.45,
	\"volumeFracMin\": 0,
	\"youngModulusMax\": 1.0,
	\"youngModulusMin\": 1e-09,
	\"materialDensity\": 0.0,
	\"a\":[[0,1]],
	\"c\":[[0.5,0.5]]
}" > tests_folder/params.json
  python main.py tests_folder/params.json tests_folder/${outputname} >output.txt
done
