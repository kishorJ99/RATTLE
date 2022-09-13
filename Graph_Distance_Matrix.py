from math import ceil, sqrt
import networkx as nx
import numpy as np
import string
import statistics

# with open("Distance Matrix.csv") as f:
#     lines = f.readlines()
#     lines = lines[1:]

#     indexs = set()

#     for i in lines:
#         line = i.split(",")
#         # print(int(line[0].split(":")[-1].strip()))
#         indexs.add(int(line[0].split(":")[-1].strip()))
#         indexs.add(int(line[1].split(":")[-1].strip()))

#     for i in range(0,146):
#         if i not in indexs:
#             print(i, " Not Found")

DM = []
NO_OF_SAMPLES = 4
DM_by_samples = []
samples_to_index = []

with open("Distance Matrix.csv") as f:
    lines = f.readlines()
    lines = lines[1:]

    indexs = set()

    noOfPoints = ceil(sqrt(len(lines)))
    
    DM = [[0 for _ in range(0, noOfPoints)] for _ in range(0, noOfPoints)]
    sample = ["" for _ in range(0, noOfPoints)]

    for i in lines:
        line = i.split(",")

        s1 = (line[0].split(":")[0].strip())
        s2 = (line[1].split(":")[0].strip())
        id1 = int(line[0].split(":")[-1].strip())
        id2 = int(line[1].split(":")[-1].strip())
        k = int(line[3].strip())
        distance = float(line[4].strip())

        
        samples_to_index.append(s1)
        samples_to_index.append(s2)

        proprtionalDistance = (distance) * 10

        sample[id1] = str(id1) + s1 
        sample[id2] = str(id2) + s2
        DM[id1][id2] = proprtionalDistance
        DM[id2][id1] = proprtionalDistance
        
        # print(int(line[0].split(":")[-1].strip()))
        # indexs.add(int(line[0].split(":")[-1].strip()))
        # indexs.add(int(line[1].split(":")[-1].strip()))

    samples_to_index = list(set(samples_to_index))

    DM_by_samples = [[[] for _ in range(0, len(samples_to_index))] for _ in range(0, len(samples_to_index))]

    for i in lines:
        line = i.split(",")

        s1 = (line[0].split(":")[0].strip())
        s2 = (line[1].split(":")[0].strip())
        id1 = int(line[0].split(":")[-1].strip())
        id2 = int(line[1].split(":")[-1].strip())
        k = int(line[3].strip())
        distance = float(line[4].strip())

        proprtionalDistance = (distance) * 10

        s1_index = samples_to_index.index(s1)
        s2_index = samples_to_index.index(s2)

        DM_by_samples[s1_index][s2_index].append(proprtionalDistance)
        DM_by_samples[s2_index][s1_index].append(proprtionalDistance)

#  Graph
dt = [('len', float)]
A = np.array(DM)
A = A.view(dt)


G = nx.from_numpy_matrix(A)
print(len(G.nodes()))
G = nx.relabel_nodes(G, dict(zip(range(len(G.nodes())), sample)))    

G = nx.drawing.nx_agraph.to_agraph(G)

G.node_attr.update(color="red", style="filled")
G.edge_attr.update(color="transparent", width="0.0")

G.draw('output/PointsGraph.png', format='png', prog='neato')

# Distance Stats
with open('output/Stats.csv', 'w') as f:
    for i in range(0, len(samples_to_index)):
        for j in range(0, i):
            line = "From: " + samples_to_index[i] + " To " + samples_to_index[j] + " Median : " + str(statistics.median(DM_by_samples[i][j])) + " Mean : " + str(statistics.mean(DM_by_samples[i][j])) + " Spread " + str(max(DM_by_samples[i][j]) - min(DM_by_samples[i][j])) + "\n"
            f.write(line)

    for i in range(0, len(samples_to_index)):
        line = "From: " + samples_to_index[i] + " To " + samples_to_index[i] + " Median : " + str(statistics.median(DM_by_samples[i][i])) + " Mean : " + str(statistics.mean(DM_by_samples[i][i])) + " Spread " + str(max(DM_by_samples[i][i]) - min(DM_by_samples[i][i])) + "\n"
        f.write(line)

