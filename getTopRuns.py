import csv
import os
from operator import itemgetter

input = "output/kmerTune_2"
output ="output/kmerTune_2/bestResults.csv"

directorys = [x[0] for x in os.walk(input) if x[0].split('/')[-1] == 'graphs']

# print (directorys)

scores =[]

for path in directorys:
    with open(path + '/../cluster_summary.tsv', 'r') as f:
        
        lines = f.readlines()
        noOfReads = len(lines)
        noOfClusters = max([int(x.split(',')[-1]) for x in lines], default=0)
        if noOfReads > 50 and noOfClusters > 2:
            print(noOfReads, noOfClusters)

            with open(path + "/ClusteringScores.txt", 'r')as f2:
                
                fileContents = f2.readline().split(',')

                # print(fileContents)

                homogeneity, completeness = float(fileContents[0].split(':')[-1]), float(fileContents[1].split(':')[-1])
                score = (path, homogeneity, completeness, noOfClusters)

                if homogeneity > 0.8 and completeness > 0.8:
                    scores.append(score)
                    # print(homogeneity, completeness)

scores.sort(key=itemgetter(1,2), reverse=True)

with open(output, 'w') as f:
    csv_out=csv.writer(f)
    csv_out.writerow(['loc','homogeneity', 'completeness', 'NO. of clusters'])
    csv_out.writerows(scores)

    # f.write('\n'.join('%s %s' % x for x in scores))