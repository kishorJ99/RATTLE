import json
import os
from Bio import SeqIO
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import completeness_score
from sklearn.metrics.cluster import homogeneity_score

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
# from subprocess import PIPE, run

# Config
configf = open('bvlenRun.json')
config = json.load(configf)

# Output
outputLocation = config['output']
os.makedirs(outputLocation, exist_ok=True)

def out(command):
    print("---------------")
    print(command)
    res = os.popen(command).read()

def addToCompositonMatrix(read, cluster, compositonMatrix,inputLabels):
    compositonMatrix[cluster][inputLabels.index(read)] += 1

def generateGraphs(location, configItem):

    # Create graphs folders
    os.makedirs(location + configItem['name'] + "/graphs/", exist_ok=True)

    # Input files - Histograms
    inputFiles = configItem['datasets']

    allLengthLists = []

    i=0

    for fileName in inputFiles:
        with open(fileName, 'r') as f:
            records = SeqIO.parse(f, 'fasta' if fileName.endswith('.fa') or fileName.endswith('.fasta') else 'fastq')
            lengthList = []

            for rec in records:
                lengthList.append(len(rec.seq))

            filePlot = sns.displot(lengthList)
            filePlot.set(title = configItem['labels'][i])
            plt.subplots_adjust(top=0.85)
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            filePlot.figure.savefig(location + configItem['name'] + '/graphs/' + configItem['labels'][i] + '_Length_Histogram.png', bbox_inches="tight")
            plt.close()

            allLengthLists += lengthList
            
            i += 1

    filePlot = sns.displot(allLengthLists)
    filePlot.set(title = 'All Lengths')
    plt.subplots_adjust(top=0.95)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    filePlot.figure.savefig(location + configItem['name'] + '/graphs/All_Length_Histogram.png', bbox_inches="tight")
    plt.close()

    # Create compositon matrix
    if len(configItem['labels']) != 0:
        clusterSummary = pd.read_csv(location + configItem['name'] + '/cluster_summary.tsv', header=None, names=['Read', 'Cluster'])

        samples = configItem['labels']

        compositonMatrix = [[ 0 for j in range (0, len(samples))] for i in range(0, max(clusterSummary['Cluster'], default=0) +1)]

        for row in zip(clusterSummary['Read'], clusterSummary['Cluster']):
            addToCompositonMatrix(row[0].split('.')[-1], row[1], compositonMatrix, configItem['labels'])

        compositonMatrixPerSample = pd.DataFrame(dict([(k, [j[i] for j in compositonMatrix]) for i,k in enumerate(samples)]), index=['Cluster ' + str(i) for i in range(0, max(clusterSummary['Cluster'], default=0) +1)])

        compositonPlt = compositonMatrixPerSample.plot(kind='bar', stacked=True)
        plt.xlabel('Clusters')
        plt.ylabel('Input Sources')
        plt.title('Cluster Compositon')
        plt.subplots_adjust(top=0.85)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        compositonPlt.figure.savefig(location + configItem['name'] + '/graphs/ClusterCompositionPlot.png', bbox_inches="tight")
        plt.close()


        # Clustering metrics
        readGroundTruth = [ configItem['labels'].index(i.split('.')[-1]) for i in clusterSummary['Read']]
        readCluster = clusterSummary['Cluster']

        homogeneity = homogeneity_score(readGroundTruth, readCluster)
        compleateness = completeness_score(readGroundTruth, readCluster)

        with open(location + configItem['name'] + '/graphs/ClusteringScores.txt', 'w') as f:
            f.write('Homogeneity : ' + str(homogeneity) + ' , Completeness : ' + str(compleateness))


# Run Rattle
for i in config['runs']:
    
    if (config['multipleOutputLoc']):
        outputLocation = config['output'] + i['output']
        os.makedirs(outputLocation, exist_ok=True)


    debug = []

    # Create cluster folders
    os.makedirs(outputLocation + i['name'] + "/clusters/", exist_ok=True)

    inputFiles = i['datasets']
    labels = i['labels']
    outputLocationForRun = outputLocation + i['name'] + '/'
    params = i['commandLineParams']


    resClusterCommand = './rattle cluster -i ' + ','.join(inputFiles) + ' -l ' + ','.join(labels) + ' -t 24 -o ' + outputLocationForRun + ' ' + ' '.join(params)
    out(resClusterCommand)

    resClusterSummaryCommand = './rattle cluster_summary -i ' + ','.join(inputFiles) + ' -l ' + ','.join(labels) + ' -c ./' + outputLocationForRun + 'clusters.out > ./' + outputLocationForRun + 'cluster_summary.tsv'
    out(resClusterSummaryCommand)

    resExtractClusterCommand = './rattle extract_clusters -i ' + ','.join(inputFiles) + ' -l ' + ','.join(labels) + ' -c ./' + outputLocationForRun + 'clusters.out -o ./' + outputLocationForRun + 'clusters'
    out(resExtractClusterCommand)

    generateGraphs(outputLocation, i)

