import json


configOutputLoc = "kmer_batch_config.json"

iso_kmer_size_start, iso_kmer_size_end = 6, 16
db_scan_mp_start, db_scan_mp_end = 3, 3
db_scan_eps_start, db_scan_eps_end = 1, 20

jsonFile = {}

jsonFile["multipleOutputLoc"] = True
jsonFile["output"] = "output/kmerTune/"
jsonFile["runs"] = []


def makeRun(name, output, datasets, labels, commandLineParams):
    run = {}
    run["name"] = name
    run["output"] = output
    run["datasets"] = datasets
    run["labels"] = labels
    run["commandLineParams"] = commandLineParams

    return run

for kmer_size in [14, 11, 16, 10, 9, 13, 15, 12, 8]:
    for eps in [12, 8, 15, 2, 17, 4, 1, 0.2, 3, 5, 10, 19, 6, 7, 9, 11, 13, 14, 16, 18]:
        
        # Ideal set
        name = "simple_bv_k" + str(kmer_size) + "_eps"+ str(eps)
        output = str(kmer_size) + "mer/"
        datasets = ["samples/sample1/50.fasta", "samples/sample1/25.fasta"]
        labels = ["S1", "S2"]
        commandLineParams = ["--iso-bitvec --iso-kmer-size " + str(kmer_size) + " --dbscan-eps " + str(eps) + " --dbscan-mp 3"]
        
        run = makeRun(name, output, datasets, labels, commandLineParams)

        jsonFile["runs"].append(run)

        # Real set
        name = "realistic_bv_k" + str(kmer_size) + "_eps"+ str(eps)
        output = str(kmer_size) + "mer/"
        datasets = ["samples/human_samples/cluster_split/468_gene:ENSG00000101335.10/0.fasta", 
                "samples/human_samples/cluster_split/468_gene:ENSG00000101335.10/1.fasta", 
                "samples/human_samples/cluster_split/317_gene:ENSG00000124172.10/0.fasta", 
                "samples/human_samples/cluster_split/274_gene:ENSG00000171858.18/0.fasta", 
                "samples/human_samples/cluster_split/269_gene:ENSG00000101210.13/0.fasta", 
                "samples/human_samples/cluster_split/269_gene:ENSG00000101210.13/1.fasta", 
                "samples/human_samples/cluster_split/269_gene:ENSG00000101210.13/2.fasta", 
                "samples/human_samples/cluster_split/253_gene:ENSG00000125868.16/0.fasta", 
                "samples/human_samples/cluster_split/253_gene:ENSG00000125868.16/1.fasta", 
                "samples/human_samples/cluster_split/125_gene:ENSG00000101439.9/0.fasta",
                "samples/human_samples/cluster_split/125_gene:ENSG00000101439.9/1.fasta",
                "samples/human_samples/cluster_split/125_gene:ENSG00000101439.9/2.fasta",
                "samples/human_samples/cluster_split/121_gene:ENSG00000101182.15/0.fasta",
                "samples/human_samples/cluster_split/121_gene:ENSG00000101182.15/1.fasta",
                "samples/human_samples/cluster_split/121_gene:ENSG00000101182.15/2.fasta",
                "samples/human_samples/cluster_split/107_gene:ENSG00000198959.12/0.fasta",
                "samples/human_samples/cluster_split/107_gene:ENSG00000198959.12/1.fasta"
            ]
        labels = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17"]
        commandLineParams = ["--iso-bitvec --iso-kmer-size " + str(kmer_size) + " --dbscan-eps " + str(eps) + " --dbscan-mp 3"]
        
        run = makeRun(name, output, datasets, labels, commandLineParams)

        jsonFile["runs"].append(run)

jsonContent = json.JSONEncoder().encode(jsonFile)

with open(configOutputLoc, "w") as f:
    f.write(jsonContent)