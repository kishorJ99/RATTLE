import json


configOutputLoc = "kmer_tune_realistic_selective_batch_config.json"

iso_kmer_size_start, iso_kmer_size_end = 4, 5
db_scan_mp_start, db_scan_mp_end = 3, 3
db_scan_eps_start, db_scan_eps_end = 1, 20

jsonFile = {}

jsonFile["multipleOutputLoc"] = True
jsonFile["output"] = "output/kmerTune_2/"
jsonFile["runs"] = []


def makeRun(name, output, datasets, labels, commandLineParams):
    run = {}
    run["name"] = name
    run["output"] = output
    run["datasets"] = datasets
    run["labels"] = labels
    run["commandLineParams"] = commandLineParams

    return run

for kmer_size in [8, 10, 13, 4, 6, 7, 14, 15, 12, 9, 16, 11, 5]:
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
        datasets = [
            "samples/human_samples/cluster_split/32_gene:ENSG00000171867.18/0.fasta",
            "samples/human_samples/cluster_split/32_gene:ENSG00000171867.18/1.fasta",
            "samples/human_samples/cluster_split/32_gene:ENSG00000171867.18/2.fasta",
            "samples/human_samples/cluster_split/32_gene:ENSG00000171867.18/3.fasta",
            ]
        labels = ["S1", "S2", "S3", "S4"]
        commandLineParams = ["--iso-bitvec --iso-kmer-size " + str(kmer_size) + " --dbscan-eps " + str(eps) + " --dbscan-mp 3"]
        
        run = makeRun(name, output, datasets, labels, commandLineParams)

        jsonFile["runs"].append(run)

jsonContent = json.JSONEncoder().encode(jsonFile)

with open(configOutputLoc, "w") as f:
    f.write(jsonContent)