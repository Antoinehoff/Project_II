#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Optimize the problem described by the given json file.
Save the resulting image at the given location.
"""

# System libs
import os
import json
import argparse
import importlib

# Third party libs
import numpy
import scipy.misc
import scipy

# Local libs
from parameters import Parameters
import multires

#Timer
import time

#plot
import matplotlib.pyplot as plt

def run(params):
    print("ndes: " + str(params.nelx) + " x " + str(params.nely))
    print("volfrac: " + str(params.volumeFracMax) + ", rmin: " +
          str(params.filterRadius) + ", penal: " + str(params.penalty))
    print("Filter method: " + params.filterType.name)

    # Load problem boundary conditions dynamically
    imported_problem = importlib.import_module("problems." + params.problemModule)
    bc = imported_problem.BoundaryConditions()

    # Solve the multires problem
    (x, nelx, nely, results, aggregated_hist) = multires.multires(params.nelx, params.nely, params, bc)

    # Convert to RGB
    x_rgb = 1.0 - numpy.reshape(numpy.repeat(x.reshape((nelx, nely)), 3, axis=1), (nelx, nely, 3))
    return (x_rgb, nelx, nely, results, aggregated_hist)

### Antoine Hoffmann 2018 EPFL
def plot_history(args,outputname):
    head, tail = os.path.split(os.path.splitext(args.input_json)[0])
    inputname ='output/'+tail+"_histories.json"
    outputname +='_histories.png'
    print('Plotting : '+ inputname +' in ' + outputname)
    aggregated_hist=[]
    with open(inputname) as json_file:
         aggregated_hist=json.load(json_file)
    y1=[]
    y2=[]
    x1=[]
    x2=[]
    temp1 = 0
    temp2 = 0
    for level in aggregated_hist.keys():
            y1.append(aggregated_hist[level]['Compliance'])
            x1.append(range(temp1,temp1+len(aggregated_hist[level]['Compliance'])))
            temp1 = temp1+len(aggregated_hist[level]['Compliance'])
            y2.append(aggregated_hist[level]['Appearance'])
            x2.append(range(temp2,temp2+len(aggregated_hist[level]['Appearance'])))
            temp2 = temp2+len(aggregated_hist[level]['Appearance'])
    plt.subplot(2,1,1)
    for i in range(len(y1)):
        plt.plot(x1[i],y1[i],label='Level '+str(i))
    plt.title('Comparison between compliance and appearance evolution')
    plt.xlabel('Step')
    plt.ylabel('Compliance Only')
    plt.legend()
    plt.subplot(2,1,2)
    for i in range(len(y2)):
        plt.plot(x2[i],y2[i],label='Level '+str(i))
    plt.xlabel('Step')
    plt.ylabel('Appearance')
    plt.legend()
    plt.savefig(outputname)

def output_name(params):
    name = 'output/'
    if params.problemType == "ProblemType.ComplianceWithSymmetry" \
                             or "ProblemType.AppearanceWithMaxComplianceAndSymmetry":
        name += 'symmetric_'
    name += params.problemModule+'_'
    name += 'mss'+str(params.maxSolverStep)+'_'
    name += str(params.nelx)+'x'+str(params.nely)+'_'
    name += 'eds'+str(params.exemplarDownsampling)
    name += '.png'
    return name
###

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_json", help="parameter file in .json format")
    return parser.parse_args()

def main(args):
    input_params = Parameters.loads(args.input_json)
    start = time.time()
    # print(os.path.splitext(args.input_json)[0] + "_results.json")
    (x_rgb_out, nelx_out, nely_out, results_out, aggregated_hist) = run(input_params)
    print(results_out, nelx_out, nely_out)
    history_name = 'output/'+os.path.split(os.path.splitext(args.input_json)[0])[1]+"_histories.json"
    with open(history_name, 'w') as json_file:
         json_file.write(json.dumps(aggregated_hist))
    outputname = output_name(input_params)
    scipy.misc.toimage(x_rgb_out.T, cmin=0.0, cmax=1.0).save(outputname+'.png')
    end = time.time()
    print('Ellapsed time : ' + str(int((end-start))/60)+'min '+str(int((end-start))%60)+'sec')

    #Analysis of the convergence:
    if input_params.record_histories == True:
        plot_history(args,outputname)

if __name__ == "__main__":
    main(parse_args())
