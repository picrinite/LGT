from JB111015.models import *

#Read.objects.all().delete()
#ReadPair.objects.all().delete()
#BlastHit.objects.all().delete()
TophatHit.objects.all().delete()

'''
# Reads table
file = open(r'DatabaseFastq.txt')

rows = file.readlines()

for i in range(0, len(rows)):
    line = rows[i]
    rowfields = line.split('\t')
    if len(rowfields) < 1:
        continue
    read = Read(
        read_id=rowfields[0],
        sample_id=Metadata.objects.get(sample_id=rowfields[1]),
        instrument_id=rowfields[2],
        run_id=rowfields[3],
        flowcell_id=rowfields[4],
        lane=rowfields[5],
        tile=rowfields[6],
        xy_coord=rowfields[7],
        pair=rowfields[8],
        seq=rowfields[9],
        seq_length=rowfields[10],
        seq_quality=rowfields[11],
        seq_category=rowfields[12]
        )
    read.save()
    if(i % 2 == 0):
        read_id_1 = rowfields[0]
    else:
        pos = rowfields[0].index('_')
        read_pair = ReadPair(
            pair_id=rowfields[0][0:pos],
            read_1=Read.objects.get(read_id=read_id_1),
            read_2=Read.objects.get(read_id=rowfields[0])
        )
        read_pair.save()

# BlastHit table
file = open(r'DatabasePairedBlast.txt')

rows = file.readlines()

for i in range(1, len(rows)):  # first line is title
    line = rows[i]
    rowfields = line.split('\t')
    if len(rowfields) < 1:
        continue
    blast_hit = BlastHit(
        read_id=Read.objects.get(read_id=rowfields[0]),
        gi_no=rowfields[1],
        accession_no=rowfields[2],
        taxon_id=rowfields[3],
        identity=rowfields[4],
        align_length=rowfields[5],
        mismatches=rowfields[6],
        evalue=rowfields[7],
        bitscore=rowfields[8],
        subject_title=rowfields[9],
        s_start=rowfields[10],
        s_end=rowfields[11],
        s_seq=rowfields[12]
        )
    blast_hit.save()

'''
file = open(r'DatabasePairedHostTophat.txt')

rows = file.readlines()

for i in range(1, len(rows)):  # first line is title
    line = rows[i]
    rowfields = line.split('\t')
    if len(rowfields) < 1:
        continue
    tophat_hit = TophatHit(
        read_id=Read.objects.get(read_id=rowfields[0]),
        sum_flags=rowfields[1],
        ref_chromosome=rowfields[2],
        chromosome_num=rowfields[3],
        mapping_position=rowfields[4],
        mapping_quality=rowfields[5],
        cigar=rowfields[6],
        overlap_status=rowfields[7],
        gene_name=rowfields[8],
        gene_id=rowfields[9]
        )
    tophat_hit.save()

