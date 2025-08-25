from sys import argv
import os
import upsetplot
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
from collections import Counter
from Bio import SeqIO
from logger import logger

def parse_representatives(filtering_dir, annotation_fn):

    repr_fn = os.path.split(annotation_fn)[1].rsplit(".",1)[0] + ".annotated.repr.tsv"
    repr_path= os.path.join(filtering_dir,repr_fn)

    representatives_set = set()
    title_line = "Cluster\tContig\tLength\tGenomad_n_genes\tGenomad_virus_score\tGenomad_taxonomy\tPalmdb_hit(bitscore >=100, pid >=50%)\tPalmdb_taxonomy\tPlant_palmprintsDB((bitscore >=100, pid >=50%))\tPalmScan positive\tEmbryophyta positive\tTop5_blastx_nr_sskingdoms\tTop5_blastx_nr_skingdoms\tBest_blastx_nr_acc\tBest_blastx_nr_taxonomy\tBest_blastx_nr_taxid\tBest_blastx_nr_title\tBest_blastx_nr_bitscore|pid\thmm|eval|bs|norm_bs|prof_cov|cont_cov\n"

    with open(annotation_fn,'r') as f, open(repr_path,'w') as out_handle:
        for line in f:
            if line.startswith('Cluster'):
                continue
            line = line.strip().split('\t')
            representatives_set.add(line[0])
        out_handle.write(title_line)
        f.seek(0)
        for line_ in f:
            if line_.startswith('Cluster'):
                continue
            new_line = line_.strip().split('\t')
            if new_line[1] in representatives_set:
                out_handle.write(line_)

    return repr_path


def remove_controls(filtering_dir, annotation_fn):

    controls_fn = os.path.split(annotation_fn)[1].rsplit('.',1)[0] + '.controls.tsv'
    controls_path = os.path.join(filtering_dir, controls_fn)

    non_controls_fn = os.path.split(annotation_fn)[1].rsplit('.',1)[0] + '.no_controls.tsv'
    non_controls_path = os.path.join(filtering_dir, non_controls_fn)

    counter = 0
    with open(annotation_fn, 'r') as f, open(non_controls_path,'w') as out_handle, open(controls_path,'w') as controls_handle:
        for line in f:
            if line.startswith('Cluster'):
                out_handle.write(line)
                controls_handle.write(line)
                continue
            new_line = line.strip().split('\t')
            if new_line[1].startswith('Neg') or new_line[1].startswith('Pos') or new_line[0].startswith('Neg') or new_line[0].startswith('Pos'):
                controls_handle.write(line)
                counter += 1
                continue
            else:
                out_handle.write(line)
    logger.log(f"Filtered out {counter} control contigs, written to {controls_path}")

    return non_controls_path


def parse_annotation_fn(annotation_fn):

    annotation_dict = {}
    with open(annotation_fn, 'r') as f:
        for line in f:
            if line.startswith('Cluster'):
                continue

            line = line.strip().split('\t')
            if line[0] not in annotation_dict:
                annotation_dict[line[0]] = [line[1:]]
            else:
                annotation_dict[line[0]].append(line[1:])

    return annotation_dict


def write_filtered_fasta(filter_outdir, annotation_fn, fasta_fn, filtered_contigs_fn= "master_contigs.filtered.fasta"):

    filtered_contigs_path =os.path.join(filter_outdir, filtered_contigs_fn)
    annotation_dict = parse_annotation_fn(annotation_fn)
    record_counter = 0
    contig_set = set()
    if not os.path.exists(filtered_contigs_path):
        for cluster, annotations in annotation_dict.items():
            for annotation in annotations:
                contig_set.add(annotation[0])
        with open(fasta_fn, 'r') as f, open(filtered_contigs_path, 'w') as out_handle:
            for record in SeqIO.parse(f, 'fasta'):
                if record.id in contig_set:
                    record_counter += 1

                    SeqIO.write(record, out_handle, 'fasta')
    logger.log(f"Written {record_counter} contigs")
    return filtered_contigs_path

def write_total_votu_contigs(filter_outdir, annotation_fn, fasta_fn, mmseqs_summary_path, filtered_contigs_fn= "master_contigs.filtered.fasta"):

    votu_contigs_path =os.path.join(filter_outdir, filtered_contigs_fn)
    annotation_dict = parse_annotation_fn(annotation_fn)
    record_counter = 0

    votu_set = set()
    contig_set = set()

    if not os.path.exists(votu_contigs_path):
        for cluster, annotations in annotation_dict.items():
            for annotation in annotations:
                votu_set.add(annotation[0])

        with open(mmseqs_summary_path, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                if line[0] in votu_set:
                    contig_set.add(line[1])

        with open(fasta_fn, 'r') as f, open(votu_contigs_path, 'w') as out_handle:
            for record in SeqIO.parse(f, 'fasta'):
                if record.id in contig_set:
                    record_counter += 1

                    SeqIO.write(record, out_handle, 'fasta')
    logger.log(f"Written {record_counter} contigs to {votu_contigs_path}")

    return votu_contigs_path







def run_mmseqs_cluster(mmseqs_inpath, mmseqs_outdir, seq_id=0.9, cov= 0.75):
    """
    Runs mmseqs easy-cluster on a given contig file and returns the path to the output file.

    :param mmseqs_inpath: Path to the contig file
    :type mmseqs_inpath: str
    :param mmseqs_outdir: Path to the output directory
    :type mmseqs_outdir: str
    :param seq_id: Sequence identity threshold, defaults to 0.9
    :type seq_id: float, optional
    :return: mmseqs_outpath: Path to the output file
    :rtype: str

    TODO: Parameters are now hardcoded, add option to change them
    """

    contig_fn = os.path.split(mmseqs_inpath)[1].split('.')[0]

    mmseqs_out_fn = f"{contig_fn}_cluster_{seq_id}.fasta"
    mmseqs_outpath = os.path.join(mmseqs_outdir, mmseqs_out_fn)

    if not os.path.exists(mmseqs_outpath):
        try:
            cmd_mmseqs = f"mmseqs easy-cluster {mmseqs_inpath}" \
                         f" {mmseqs_outpath} {os.path.join(mmseqs_outdir, 'tmp')} " \
                         f"--min-seq-id {seq_id} -c {cov} --cov-mode 1 --cluster-mode 2"

            subprocess.check_output(cmd_mmseqs,
                                    shell=True)

        except(subprocess.CalledProcessError):
            logger.log(f"Something went wrong with mmseqs for sample {mmseqs_out_fn}")
           # logger.log(
           #      f"Something went wrong with mmseqs for sample {mmseqs_out_fn}")

    for file in os.listdir(mmseqs_outdir):
        if file.endswith('cluster.tsv'):
            mmseqs_summary = file
            # logger.log(f"mmseqs cluster file found: {file}")
        else:
            continue
            #logger.log(f"Not the cluster file I am looking for: {file}")

    mmseqs_summary_path = os.path.join(mmseqs_outdir, mmseqs_summary)

    return mmseqs_summary_path


def combine_cluster_annotations(mmseqs_summary_path, annotation_fn, combined_fn):

    annot_dict_clust = parse_annotation_fn(annotation_fn)
    annot_dict_cont = {}

    for cluster, annotations in annot_dict_clust.items():
        for annotation in annotations:
            if annotation[0] not in annot_dict_cont:
                annot_dict_cont[annotation[0]] = [annotation[1:]]
            else:
                annot_dict_cont[annotation[0]].append(annotation[1:])

    with open(mmseqs_summary_path, 'r') as f, open(combined_fn, 'w') as out_handle:
        for line in f:
            line = line.strip().split('\t')
            cluster = line[0]
            contig = line[1]
            if contig in annot_dict_cont.keys():
                for annotation in annot_dict_cont[contig]:
                    annot = '\t'.join(annotation)
                    out_handle.write(f"{cluster}\t{contig}\t{annot}\n")
            else:
                continue

    return combined_fn


def cluster_stats(combined_fn):

    combined_dict = parse_annotation_fn(combined_fn)

    cluster_stats = {}
    for cluster, annotations in combined_dict.items():
        cluster_stats[cluster] = len(annotations)

    logger.log(f"Total number of filtered clusters: {len(cluster_stats)} with 0.9 seq_id and 0.75 cov cutoff")
    with open("cluster_size.tsv", 'w') as out_handle:
        for cluster, size in cluster_stats.items():
            out_handle.write(f"{cluster}\t{size}\n")

    # Create histogram
    cluster_stats = pd.DataFrame.from_dict(cluster_stats, orient='index')
    cluster_stats = cluster_stats.reset_index()
    cluster_stats.columns = ['Cluster', 'Number of contigs in cluster']
    plt.hist(cluster_stats['Number of contigs in cluster'], bins=100)
    plt.xlabel('Number of contigs in cluster')
    plt.ylabel('Number of clusters')
    plt.title('Cluster size distribution')
    plt.savefig('cluster_size_distribution.png')
    return cluster_stats


def cluster_filtered(filter_outdir,annotation_fn, mmseqs_inpath, mmseqs_outdir, combined_fn= "master_contigs.annotated.no_eves.clustered.tsv"):

    combined_fn_path = os.path.join(filter_outdir, combined_fn)
    filtered_fasta = write_filtered_fasta(filter_outdir,annotation_fn, mmseqs_inpath)
    mmseqs_summary_path = run_mmseqs_cluster(filtered_fasta, mmseqs_outdir)
    combined_fn_path = combine_cluster_annotations(mmseqs_summary_path, annotation_fn, combined_fn_path)

    return combined_fn_path, mmseqs_summary_path


def filter_embryophyta(annotation_dict):

    embryophyta_dict = {}
    filtered_embryophyta_dict = {}
    for cluster, annotations in annotation_dict.items():
        for annotation in annotations:
            if annotation[9] != 'NA':
                annot_line = annotation[9].split('|')
                for annot in annot_line:
                    if annot.startswith('bs:'):
                        bitscore = float(annot.split(':')[1])
                    elif annot.startswith('eval:'):
                        evalue = float(annot.split(':')[1])
                    elif annot.startswith('pid:'):
                        pid = float(annot.split(':')[1])
                    elif annot.startswith('query_cov:'):
                        query_cov = float(annot.split(':')[1])
                if pid >= 50 and query_cov >= 0.25:
                    if cluster not in embryophyta_dict:
                        embryophyta_dict[cluster] = [annotation]
                    else:
                        embryophyta_dict[cluster].append(annotation)
                else:
                    if cluster not in filtered_embryophyta_dict:
                        filtered_embryophyta_dict[cluster] = [annotation]
                    else:
                        filtered_embryophyta_dict[cluster].append(annotation)
            else:
                if cluster not in filtered_embryophyta_dict:
                    filtered_embryophyta_dict[cluster] = [annotation]
                else:
                    filtered_embryophyta_dict[cluster].append(annotation)


    return embryophyta_dict, filtered_embryophyta_dict


def prefilter_mmseqs_nt_non_virus(mmseqs_nt_fn, mmseqs_nt_filtered_fn, coverage = 0.25, seq_id_thres= 50):

    with open(mmseqs_nt_fn, 'r') as f, open(mmseqs_nt_filtered_fn, 'w') as out_handle:
        for line in f:
            line = line.strip().split('\t')
            lineage = line[16].split(';')
            lineage_tag= lineage[0]
            if float(line[12]) >= coverage and float(line[2]) >= seq_id_thres and lineage_tag == '-_cellular organisms':
                out_handle.write('\t'.join(line)+'\n')
    return mmseqs_nt_filtered_fn

def prefilter_blastx_nr_non_virus(blastx_nr_fn, blastx_nr_filtered_fn, coverage = 25, seq_id_thres= 50):

    with open(blastx_nr_fn, 'r') as f, open(blastx_nr_filtered_fn, 'w') as out_handle:
        for line in f:
            line = line.strip().split('\t')
            if float(line[18]) >= coverage and float(line[2]) >= seq_id_thres and 'viruses' not in line[14].lower():
                out_handle.write('\t'.join(line)+'\n')
    return blastx_nr_filtered_fn

def filter_mmseqs_nt_non_virus(embr_filtered_dict,embr_dict,blastx_nr_filtered_fn):

    blastx_non_vir_set = set()

    with open(blastx_nr_filtered_fn, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            blastx_non_vir_set.add(line[0])

    virus_dict = {}

    for cluster, annotations in embr_filtered_dict.items():
        for annotation in annotations:
            if annotation[0] in blastx_non_vir_set:
                if cluster not in embr_dict:
                    embr_dict[cluster] = [annotation]
                else:
                    embr_dict[cluster].append(annotation)
            else:
                if cluster not in virus_dict:
                    virus_dict[cluster] = [annotation]
                else:
                    virus_dict[cluster].append(annotation)

    return virus_dict, embr_dict


def filter_eves(filtering_dir, annotation_fn, mmseqs_nt_out):

    eves_fn = os.path.join(filtering_dir,"eve_candidates.tsv")
    non_eves_fn = os.path.join(filtering_dir, "eve_filtered.tsv")


    annot_dict = parse_annotation_fn(annotation_fn)
    embryophyta_dict, filtered_embryophyta_dict = filter_embryophyta(annot_dict)

    # prefilter blastx_nr for non-virus
    mmseqs_nt_filtered_fn = os.path.split(mmseqs_nt_out)[1].rsplit(".",1)[0] + "mmseqs_nt_filtered_non_virus.m8"
    mmseqs_nt_filtered_path = os.path.join(filtering_dir,mmseqs_nt_filtered_fn)

    # Prefilter mmseqs output
    mmseqs_nt_filtered_path = prefilter_mmseqs_nt_non_virus(mmseqs_nt_out, mmseqs_nt_filtered_path)

    # filter non-embryophyta dict for non-virus blastx_nr contigs
    virus_dict, non_virus_dict = filter_mmseqs_nt_non_virus(filtered_embryophyta_dict,embryophyta_dict, mmseqs_nt_filtered_path)

    eve_counter = 0

    with open(eves_fn,'w') as out_handle:
        for cluster, annotations in non_virus_dict.items():
            for annotation in annotations:
                eve_counter += 1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    virus_counter = 0

    with open(non_eves_fn,'w') as out_handle:
        for cluster, annotations in virus_dict.items():
            for annotation in annotations:
                virus_counter += 1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    logger.log(f"Filtered out {eve_counter} EVE contigs")
    logger.log(f"{virus_counter} non-EVE contigs survived filtering!")

    return eves_fn, non_eves_fn



def filter_palmscan(filter_outdir, virus_out, palmscan_out= "master_cluster.annotated.no_eves.palmscan.tsv", non_palm_out= "master_cluster.annotated.no_eves.no_palmscan.tsv"):

    palmscan_dict = {}
    non_palm_dict = {}

    palmscan_out_fn = os.path.split(virus_out)[1].rsplit(".")[0] + "no_eves.palmscan.tsv"
    palmscan_out_path = os.path.join(filter_outdir, palmscan_out_fn)

    non_palmscan_out_fn = os.path.split(virus_out)[1].rsplit(".")[0] + "no_eves.no_palmscan.tsv"
    non_palmscan_out_path = os.path.join(filter_outdir,non_palmscan_out_fn)


    virus_out_dict = parse_annotation_fn(virus_out)

    for cluster, annotations in virus_out_dict.items():
        for annotation in annotations:
            if annotation[8] == 'True':

                if cluster not in palmscan_dict:
                    palmscan_dict[cluster] = [annotation]
                else:
                    palmscan_dict[cluster].append(annotation)
            else:
                if cluster not in non_palm_dict:
                    non_palm_dict[cluster] = [annotation]
                else:
                    non_palm_dict[cluster].append(annotation)

    palmscan_counter = 0
    with open(palmscan_out_path,'w') as out_handle:
        for cluster, annotations in palmscan_dict.items():
            for annotation in annotations:
                palmscan_counter += 1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    non_palm_counter = 0
    with open(non_palmscan_out_path,'w') as out_handle:
        for cluster, annotations in non_palm_dict.items():
            for annotation in annotations:
                non_palm_counter += 1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    logger.log(f"{non_palm_counter} contigs filtered due to no palmprint hits")
    logger.log(f"{palmscan_counter} contigs survived palmscan filtering")

    return palmscan_out_path, non_palmscan_out_path


def filter_prof_cov(filter_outdir, palmscan_out,non_palmscan_out,prof_cov_thresh=0.75):


    prof_cov_out_fn = os.path.split(palmscan_out)[1].rsplit(".",1)[0] + ".prof_cov.tsv"
    prof_cov_out_path = os.path.join(filter_outdir, prof_cov_out_fn)

    non_prof_cov_out_fn = os.path.split(non_palmscan_out)[1].rsplit(".",1)[0] + ".non_prof_cov.tsv"
    non_prof_cov_out_path = os.path.join(filter_outdir, non_prof_cov_out_fn)


    non_palm_dict = parse_annotation_fn(non_palmscan_out)
    palmscan_dict = parse_annotation_fn(palmscan_out)
    prof_cov_dict = {}

    ambig_dict = {}
    for cluster, annotations in non_palm_dict.items():
        for annotation in annotations:
            annotation_line = annotation[17].split('|')
            for annot in annotation_line:
                if annot.startswith('bitscore:'):
                    bitscore = float(annot.split(':')[1])
                elif annot.startswith('eval:'):
                    evalue = float(annot.split(':')[1])
                elif annot.startswith('norm_bitscore:'):
                    norm_bitscore = float(annot.split(':')[1])
                elif annot.startswith('profile_cov:'):
                    prof_cov = float(annot.split(':')[1])
                elif annot.startswith('contig_cov:'):
                    contig_cov = float(annot.split(':')[1])

            if prof_cov >= prof_cov_thresh:
                if cluster not in palmscan_dict:
                    palmscan_dict[cluster] = [annotation]
                else:
                    palmscan_dict[cluster].append(annotation)
            else:
                if cluster not in ambig_dict:
                    ambig_dict[cluster] = [annotation]
                else:
                    ambig_dict[cluster].append(annotation)

    prof_counter = 0
    with open(prof_cov_out_path,'w') as out_handle:
        for cluster, annotations in palmscan_dict.items():
            for annotation in annotations:
                prof_counter +=1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    non_prof_counter = 0
    with open(non_prof_cov_out_path,'w') as out_handle:
        for cluster, annotations in ambig_dict.items():
            for annotation in annotations:
                non_prof_counter +=1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    logger.log(f"{non_prof_counter} contigs filtered due to hmm profile coverage less than {prof_cov_thresh}%")
    logger.log(f"{prof_counter} contigs survived hmm filtering ")

    return prof_cov_out_path, non_prof_cov_out_path


def filter_palmdb(filter_outdir, prof_cov_out, non_prof_cov_out, palmdb_thresh=90):

    palmdb_out_fn = os.path.split(prof_cov_out)[1].rsplit(".",1)[0] + ".palmdb.tsv"
    palmdb_out_path = os.path.join(filter_outdir, palmdb_out_fn)

    non_palmdb_out_fn = os.path.split(non_prof_cov_out)[1].rsplit(".",1)[0] + ".non_palmdb.tsv"
    non_palmdb_out_path = os.path.join(filter_outdir, non_palmdb_out_fn)

    prof_cov_dict = parse_annotation_fn(prof_cov_out)
    non_prof_cov_dict = parse_annotation_fn(non_prof_cov_out)
    non_palmdb_dict = {}

    for cluster, annotations in non_prof_cov_dict.items():
        for annotation in annotations:
            if annotation[5] != 'NA':
                annot_line = annotation[5].split('|')
                for annot in annot_line:
                    if annot.startswith('bs:'):
                        bitscore = float(annot.split(':')[1])
                    elif annot.startswith('pid:'):
                        pid = float(annot.split(':')[1])
                    elif annot.startswith('tcov:'):
                        tcov = float(annot.split(':')[1])
                if tcov >= palmdb_thresh:

                    if cluster not in prof_cov_dict:
                        prof_cov_dict[cluster] = [annotation]
                    else:
                        prof_cov_dict[cluster].append(annotation)
                else:
                    if cluster not in non_palmdb_dict:
                        non_palmdb_dict[cluster] = [annotation]
                    else:
                        non_palmdb_dict[cluster].append(annotation)
            else:
                if cluster not in non_palmdb_dict:
                    non_palmdb_dict[cluster] = [annotation]
                else:
                    non_palmdb_dict[cluster].append(annotation)

    palmdb_counter = 0
    non_palmdb_counter = 0
    with open(palmdb_out_path,'w') as out_handle:
        for cluster, annotations in prof_cov_dict.items():
            for annotation in annotations:
                palmdb_counter += 1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    with open(non_palmdb_out_path,'w') as out_handle:
        for cluster, annotations in non_palmdb_dict.items():
            for annotation in annotations:
                non_palmdb_counter += 1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    logger.log(f"{palmdb_counter} contigs survived palmdb filtering")
    logger.log(f"{non_palmdb_counter} contigs filtered due to less than {palmdb_thresh}% coverage to palmprint ")

    return palmdb_out_path, non_palmdb_out_path


def filter_blastxnr_virus(palmdb_out, non_palmdb_out, blastx_nr_out= "master_contigs.annotated.no_eves.palmscan.prof_cov.palmdb.blastx_nr_virus.tsv", non_blastx_nr_out= "master_contigs.annotated.non_blastx_nr.tsv"):

    palmdb_dict = parse_annotation_fn(palmdb_out)
    non_palmdb_dict = parse_annotation_fn(non_palmdb_out)

    non_blastx_nr_dict = {}


    for cluster, annotations in non_palmdb_dict.items():
        for annotation in annotations:
            if annotation[16] != 'NA' and "viruses" in annotation[10].split(",")[0].lower():
                annot_line = annotation[16].split('|')
                for annot in annot_line:
                    if annot.startswith('bs:'):
                        bitscore = float(annot.split(':')[1])
                    elif annot.startswith('pid:'):
                        pid = float(annot.split(':')[1])
                    elif annot.startswith('best_hit_cov:'):
                        best_hit_cov = float(annot.split(':')[1])
                if best_hit_cov >= 80:

                    if cluster not in palmdb_dict:
                        palmdb_dict[cluster] = [annotation]
                    else:
                        palmdb_dict[cluster].append(annotation)
                else:
                    if cluster not in non_blastx_nr_dict:
                        non_blastx_nr_dict[cluster] = [annotation]
                    else:
                        non_blastx_nr_dict[cluster].append(annotation)
            else:
                if cluster not in non_blastx_nr_dict:
                    non_blastx_nr_dict[cluster] = [annotation]
                else:
                    non_blastx_nr_dict[cluster].append(annotation)

    blastx_virus_counter = 0
    non_blastx_virus_counter = 0
    with open(blastx_nr_out,'w') as out_handle:
        for cluster, annotations in palmdb_dict.items():
            for annotation in annotations:
                blastx_virus_counter += 1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    with open(non_blastx_nr_out,'w') as out_handle:
        for cluster, annotations in non_blastx_nr_dict.items():
            for annotation in annotations:
                non_blastx_virus_counter += 1
                annotation = '\t'.join(annotation)
                out_handle.write(f"{cluster}\t{annotation}\n")

    print(f"{blastx_virus_counter} survived filtering for blastx_nr_virus")
    print(f"{non_blastx_virus_counter} filtered out due to less that 80% coverage with a viral blastx hit")
    return blastx_nr_out, non_blastx_nr_out


def filter_blastx_non_virus_2(embryo_filt_dict, embryo_dict):

    virus_dict = {}

    for cluster, annotations in embryo_filt_dict.items():
        for annotation in annotations:
            if annotation[16] != 'NA':
                if annotation[10].split(",")[0] != "Viruses":
                    annot_line = annotation[16].split('|')
                    for annot in annot_line:
                        if annot.startswith('bs:'):
                            bitscore = float(annot.split(':')[1])
                        elif annot.startswith('pid:'):
                            pid = float(annot.split(':')[1])
                        elif annot.startswith('best_hit_cov:'):
                            query_cov = float(annot.split(':')[1])
                    if pid >= 0 and query_cov >= 50:
                        if cluster not in embryo_dict:
                            embryo_dict[cluster] = [annotation]
                        else:
                            embryo_dict[cluster].append(annotation)
                    else:
                        if cluster not in virus_dict:
                            virus_dict[cluster] = [annotation]
                        else:
                            virus_dict[cluster].append(annotation)
                        continue
                else:
                    if cluster not in virus_dict:
                        virus_dict[cluster] = [annotation]
                    else:
                        virus_dict[cluster].append(annotation)
                    continue
            else:
                if cluster not in virus_dict:
                    virus_dict[cluster] = [annotation]
                else:
                    virus_dict[cluster].append(annotation)
                continue

    return virus_dict, embryo_dict

def form_upset_dict(annotation_dict):

    upset_dict = {}
    upset_dict_formatted = {}
    upset_dict_length = {}

    for cluster, annotations in annotation_dict.items():
        upset_dict_formatted[cluster] = []
        tmp_formatted_list = []

        upset_dict[cluster] = {}
        upset_dict[cluster]['Genomad'] = False
        upset_dict[cluster]['Palmdb'] = False
        upset_dict[cluster]['Plant_Palmprints'] = False
        upset_dict[cluster]['PalmScan'] = False
        upset_dict[cluster]['BlastX_nr_top_5_virus'] = False

        for annotation_cols in annotations:
            genoma_vir_score = annotation_cols[3]
            palmdb_hit = annotation_cols[5]
            plant_palmprint_hit = annotation_cols[7]
            palmscan_bool = annotation_cols[8]
            sskingdom_list = annotation_cols[10].split(',')
            if genoma_vir_score != 'NA':
                upset_dict[cluster]['Genomad'] = True
                tmp_formatted_list.append('Genomad')


            if palmdb_hit != 'NA':
                upset_dict[cluster]['Palmdb'] = True
                tmp_formatted_list.append('Palmdb')

            if plant_palmprint_hit != 'NA':
                upset_dict[cluster]['Plant_Virus_Palmprints'] = True
                tmp_formatted_list.append('Plant_Virus_Palmprints')

            if palmscan_bool == 'True':
                upset_dict[cluster]['PalmScan'] = True
                tmp_formatted_list.append('PalmScan')

            if "Viruses" in sskingdom_list:
                upset_dict[cluster]['BlastX_nr_top_5_virus'] = True
                tmp_formatted_list.append('BlastX_nr_top_5_virus')

        set_formatted = set(tmp_formatted_list)
        upset_dict_formatted[cluster] = list(set_formatted)
        upset_dict_length[cluster] = len(set_formatted)

    lists = [lst for lst in upset_dict_formatted.values()]
    upset_data = upsetplot.from_memberships(lists)
    ax_dict = upsetplot.UpSet(upset_data, subset_size='count', show_counts=True, sort_by= "degree").plot()
    plt.title('Upset plot of clusters with a given number of methods of annotation')
    plt.savefig('upset_plot.png', dpi=300)
    plt.close()

    counter = Counter(upset_dict_length.values())
    bars = plt.bar(counter.keys(), counter.values())
    for bar, count in zip(bars, counter.values()):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), str(count), ha='center', va='bottom')

    plt.title('Number of clusters with a given number of methods of annotation')
    plt.xlabel('Number of methods of annotation')
    plt.ylabel('Number of clusters')
    plt.savefig('bar_plot.png', dpi=300)
    plt.close()

    return upset_dict

def create_upsetplot(annotation_fn):

    annotation_dict = parse_annotation_fn(annotation_fn)
    upset_dict = form_upset_dict(annotation_dict)
    upset_df = pd.DataFrame.from_dict(upset_dict)
    transposed_df = upset_df.transpose()
    transposed_df.to_csv('upset_df_filtered.tsv', sep='\t')

    return None


def wrap_upsetplot(annotation_fn):

        annotation_dict = parse_annotation_fn(annotation_fn)
        create_upsetplot(annotation_dict)
        return None


def wrap_filtering(overview_dir, annotation_path,mmseqs_nt_out_path,contig_path):


    filter_outdir = os.path.join(overview_dir, "filtering")
    if not os.path.exists(filter_outdir):
        os.makedirs(filter_outdir)


    embryophyta_out = "master_contigs.annotated.embryophyta.tsv"
    # Filter out control contigs
    non_controls_fn = remove_controls(filter_outdir,annotation_path)
    # Filter out eve candidate contigs
    # prefilter_mmseqs_nt_non_virus("master_contigs.mmseqs_nt.m8", "master_contigs.mmseqs_nt_filtered_non_virus.m8")
    eves_fn, non_eves_fn = filter_eves(filter_outdir, non_controls_fn, mmseqs_nt_out_path)
    # Cluster the filtered contigs
    mmseqs_outdir = os.path.join(filter_outdir, "mmseqs_cluster_filtered")
    if not os.path.exists(mmseqs_outdir):
        os.makedirs(mmseqs_outdir)

    combined_fn, mmseqs_summary_path = cluster_filtered(filter_outdir,non_eves_fn, contig_path, mmseqs_outdir)

    # Parse the combined annotation file reprmesentatives
    virus_out_repr = parse_representatives(filter_outdir,combined_fn)
    # Run the filtering steps

    palmscan_out, non_palmscan_out = filter_palmscan(filter_outdir,virus_out_repr)
    prof_cov_out, non_prof_cov_out = filter_prof_cov(filter_outdir,palmscan_out,non_palmscan_out)
    palmdb_out, non_palmdb_out = filter_palmdb(filter_outdir,prof_cov_out, non_prof_cov_out)
    filtered_contigs_path = "vOTUs.fasta"
    write_filtered_fasta(filter_outdir,palmdb_out,contig_path, filtered_contigs_path)

    votu_total_contigs_path = "vOTUs.total_contigs.fasta"
    write_total_votu_contigs(filter_outdir,palmdb_out,contig_path, mmseqs_summary_path, votu_total_contigs_path)
    # cluster_stats(palmdb_out)
    # create_upsetplot(blastx_nr_out)



if __name__ == "__main__":
    # annotation_fn = argv[1]
    # blastx_nr_fn = argv[2]
    # wrap_filtering(overview_dir, annotation_fn,mmseqs_nt_out_path, contig_path)
    pass