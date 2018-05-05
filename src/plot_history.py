#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Plot the tracking of the compliance and appearance evolution during an
optimization process
"""
import os
import json
import argparse
import numpy
import matplotlib.pyplot as plt

### Antoine Hoffmann 2018 EPFL
def plot_history(aggregated_hist, mss, filename):
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

    plt.savefig(filename)

###

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_json")
    return parser.parse_args()

def main(args):
    inputname = args.input_json
    outputname = os.path.splitext(args.input_json)[0]+'.png'
    print('Plotting : '+ inputname +' in ' + outputname)
    aggregated_hist=[]
    with open(inputname) as json_file:
         aggregated_hist=json.load(json_file)
    mss = len(aggregated_hist['level0']['Compliance'])
    plot_history(aggregated_hist,mss,outputname)
    for level in aggregated_hist.keys():
        for index in aggregated_hist[level].keys():
            list = aggregated_hist[level][index]
            print(' length of '+index+' at '+level+' = '+str(len(list)))


if __name__ == "__main__":
    main(parse_args())
