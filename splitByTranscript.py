from Bio import SeqIO
from Bio.SeqIO import FastaIO
import os

sampleFileLoc = 'samples/human_samples/sample.fa'
refFileLoc = 'samples/human_samples/ref.fa'
outputLoc = 'samples/human_samples/cluster_split/'

transcriptIds = []
seqRecords = []

with open(sampleFileLoc, 'r') as f:
    records = SeqIO.parse(f, 'fasta' if sampleFileLoc.endswith('.fa') or sampleFileLoc.endswith('.fasta') else 'fastq')

    for rec in records:
        transcriptId = rec.id.split(',')[-1]

        if transcriptId in transcriptIds:
            i = transcriptIds.index(transcriptId)
            seqRecords[i].append(rec)

        else:
            transcriptIds.append(transcriptId)
            seqRecords.append([rec])

geneIds = []
seqRecordsGroupedByGene =[]
doneList=[]

with open(refFileLoc, 'r') as f:
    refFileRecords = SeqIO.parse(f, 'fasta' if refFileLoc.endswith('.fa') or refFileLoc.endswith('.fasta') else 'fastq')

    for rec in refFileRecords:
        headerRow = rec.description.split(' ')
        print(headerRow)
        transcriptId = headerRow[0]
        geneId = headerRow[3]

        if not transcriptId in doneList:
            if transcriptId in transcriptIds:
                if geneId in geneIds:
                    seqRecordsGroupedByGene[geneIds.index(geneId)].append(seqRecords[transcriptIds.index(transcriptId)])
                else:
                    geneIds.append(geneId)
                    seqRecordsGroupedByGene.append([seqRecords[transcriptIds.index(transcriptId)]])
            doneList.append(transcriptId)

for i, geneGroup in enumerate(seqRecordsGroupedByGene):

    geneGroupFolder = outputLoc + str(len(geneGroup[0])) + '_' + geneIds[i]

    os.makedirs(geneGroupFolder, exist_ok=True)

    for j, transcriptGroup in enumerate(geneGroup):
        with open(geneGroupFolder + '/' + str(j) + ".fasta", "w") as output_handle:
            fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
            fasta_out.write_file(transcriptGroup)



# for i,v in enumerate(seqRecords):
#     with open(outputLoc + str(i) + ".fasta", "w") as output_handle:
#         fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
#         fasta_out.write_file(v)
            