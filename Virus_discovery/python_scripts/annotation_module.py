import os
import subprocess
from logger import logger



### DATABASES ###
PLANT_VIRUS_DB = "/lustre/BIF/nobackup/karap002/PHD/DBs/plant_virus_database/palmscan_plant_vir_db_palmprints.dmnd"
PALM_DB = "/lustre/BIF/nobackup/karap002/PHD/DBs/palm_db/palm_db_uniques.dmnd"
PALM_DB_TAX = "/lustre/BIF/nobackup/karap002/PHD/DBs/palm_db/u_tax.tsv"
BLAST_NR_DB = "/lustre/BIF/nobackup/karap002/PHD/DBs/blast_nr/diamond_nr.dmnd"
GENOMAD_DB = "/lustre/BIF/nobackup/karap002/PHD/DBs/genomad_db"
EMBRYOPHYTA_DB = "/lustre/BIF/nobackup/karap002/PHD/DBs/viridiplantae_nt/embryophyta_refseq_plus_65561_39414_20231107"
MMSEQS_NT= "/lustre/BIF/nobackup/karap002/PHD/DBs/mmseqs_nt/mmseqs_nt"

def run_mmseqs_cluster(mmseqs_inpath, mmseqs_outdir, seq_id=0.9, cov= 0.8):
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
           logger.log(
                f"Something went wrong with mmseqs for sample {mmseqs_out_fn}")


    for file in os.listdir(mmseqs_outdir):
        if file.endswith('cluster.tsv'):
            mmseqs_summary = file
            logger.log(f"mmseqs cluster file found: {file}")
        else:
            logger.log(f"Not the cluster file I am looking for: {file}")

    mmseqs_summary_path = os.path.join(mmseqs_outdir, mmseqs_summary)


    return mmseqs_summary_path



def run_simple_dblastx(blast_outpath, contig_path, db_path):
    """
    Runs diamond blastx on a given contig file and returns the path to the output file.

    :param blast_outpath: Path to the output file
    :type blast_outpath: str
    :param contig_path: Path to the contig file
    :type contig_path: str
    :param db_path: Path to the database
    :type db_path: str
    :return: blast_outpath: Path to the output file
    :rtype: str

    """

    if not os.path.exists(blast_outpath):

        try:
            # cmd_diamond = f"diamond blastx -q {contig_path} -d {db_path} " \
            #               f"--outfmt 6 --threads 20 --max-target-seqs 5 -e 0.001 --out {blast_outpath}"
            cmd_diamond = f"diamond blastx -q {contig_path} -d {db_path} " \
                          f"--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue" \
                          f" bitscore stitle qcovhsp scovhsp slen qlen " \
                          f"--threads 20 --max-target-seqs 5 -e 0.001 --out {blast_outpath}"

            subprocess.check_output(cmd_diamond, shell=True)
        except(subprocess.CalledProcessError):

           logger.log(
                f"Something went wrong with diamond for sample:{os.path.split(blast_outpath)[0]}")

    return blast_outpath

def run_dblasx(blast_outpath, contig_path, db_path):
    """
    Runs diamond blastx on a given contig file and returns the path to the output file.

    :param blast_outpath: Path to the output file
    :type blast_outpath: str
    :param contig_path: Path to the contig file
    :type contig_path: str
    :param db_path: Path to the database
    :type db_path: str
    :return: blast_outpath: Path to the output file
    :rtype: str

    Blastx parameters as in Dieke's pipeline, except the E-val that is set to 0.001
    """

    if not os.path.exists(blast_outpath):

        try:
            cmd_diamond = f"diamond blastx -q {contig_path} -d {db_path} " \
                          f"--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue" \
                          f" bitscore staxids sscinames sskingdoms skingdoms sphylums stitle qcovhsp slen qlen " \
                          f"--threads 20 --max-target-seqs 5 -e 0.001 --out {blast_outpath}"
            subprocess.check_output(cmd_diamond,shell=True)
        except(subprocess.CalledProcessError):

           logger.log(
                f"Something went wrong with diamond for sample:{os.path.split(blast_outpath)[0]}")


    return blast_outpath


def run_palmscan(sample, palmscan_outdir, master_contigs_path):
    """
    Runs palmscan on a given contig file and returns the path to the output file.

    :param sample: Name of the sample
    :type sample: str
    :param overview_outdir_path: Path to the overview output directory
    :type overview_outdir_path: str
    :param master_contigs_path: Path to the master contig file
    :return: fev_fn: Path to the palmscan output file with the .fev extension

    TODO: PalmScan is a hard dependency, it needs to be downloaded manually as it is not a conda package.
    TODO: As it is not crucial and is used only for validation, it may be removed in the future.
    """


    if not os.path.exists(palmscan_outdir):
        os.makedirs(palmscan_outdir)

    report_fn = f'{palmscan_outdir}/{sample}_palmscan_report.tsv'
    fev_fn = f'{palmscan_outdir}/{sample}_palmscan.fev'
    fasta_fn = f'{palmscan_outdir}/{sample}_palmscan.fasta'

    try:
        cmd_palmscan = f"/lustre/BIF/nobackup/karap002/PHD/palmscan/palmscan/bin/palmscan -search_pp " \
                       f"{master_contigs_path} -report {report_fn} -ppout {fasta_fn} -fevout {fev_fn} -rdrp -loconf"
        subprocess.check_output(cmd_palmscan, shell=True)
    except:
       logger.log(f"Something went wrong with palmscan for sample:{sample}")

    return fev_fn


def mmseqs_easy_search_simple(master_contigs,db, db_name, mmseqs_outdir):
    """    Runs mmseqs search on a given contig file and returns the path to the output file.

    :param master_contigs: Path to the master contig file
    :type master_contigs: str
    :return: mmseqs_outpath: Path to the output file
    :rtype: str
    """

    contig_fn = os.path.split(master_contigs)[1].split('.')[0]
    mmseqs_out_fn = f"{contig_fn}_{db_name}_search.m8"
    mmseqs_outpath = os.path.join(mmseqs_outdir, mmseqs_out_fn)
    if not os.path.exists(mmseqs_outpath):
        try:
            cmd_mmseqs = f"mmseqs easy-search --max-seqs 1 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov --search-type 3 {master_contigs} {db} {mmseqs_outpath}" \
                         f" {os.path.join(mmseqs_outdir, 'tmp')} "
            subprocess.check_output(cmd_mmseqs,
                                    shell=True)
        except(subprocess.CalledProcessError):
           logger.log(
                f"Something went wrong with mmseqs for sample {mmseqs_out_fn}")

    return mmseqs_outpath

def mmseqs_easy_search_tax(master_contigs,db, db_name, mmseqs_outdir):
    """    Runs mmseqs search on a given contig file and returns the path to the output file. Requires seqTaxDB.

    :param master_contigs: Path to the master contig file
    :type master_contigs: str
    :return: mmseqs_outpath: Path to the output file
    :rtype: str
    """

    contig_fn = os.path.split(master_contigs)[1].split('.')[0]
    mmseqs_out_fn = f"{contig_fn}.{db_name}.m8"
    mmseqs_outpath = os.path.join(mmseqs_outdir, mmseqs_out_fn)
    if not os.path.exists(mmseqs_outpath):
        try:
            cmd_mmseqs = f"mmseqs easy-search --max-seqs 5 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov,taxid,taxname,taxlineage -e 0.001 {master_contigs} {db} {mmseqs_outpath}" \
                         f" {os.path.join(mmseqs_outdir, 'tmp')} "
            subprocess.check_output(cmd_mmseqs,
                                    shell=True)
        except(subprocess.CalledProcessError):
           logger.log(
                f"Something went wrong with mmseqs for sample {mmseqs_out_fn}")

    return mmseqs_outpath

def run_genomad(master_contigs,genomad_db,genomad_outdir):

    if not os.path.exists(genomad_outdir):
        os.makedirs(genomad_outdir)

    genomad_out_path = os.path.join(genomad_outdir, f"{os.path.split(master_contigs)[1].split('.')[0]}_summary", f"{os.path.split(master_contigs)[1].split('.')[0]}_virus_summary.tsv")

    if not os.path.exists(genomad_out_path):
        try:
            cmd_genomad = f"conda run -n genomad genomad end-to-end --disable-nn-classification {master_contigs} {genomad_outdir} {genomad_db}"
            subprocess.check_output(cmd_genomad, shell=True)
        except:
            logger.log(f"Something went wrong with genomad for sample:{master_contigs}")


    return genomad_out_path


def wrap_annotation(master_contigs,overview_outdir):

    if not os.path.exists(os.path.join(overview_outdir,"annotate")):
        os.makedirs(os.path.join(overview_outdir,"annotate"))

    # Run mmseqs cluster
    mmseqs_outdir_clust = os.path.join(overview_outdir,"annotate","mmseqs_cluster")
    if not os.path.exists(mmseqs_outdir_clust):
        os.makedirs(mmseqs_outdir_clust)
    mmseqs_outpath_clust = run_mmseqs_cluster(master_contigs, mmseqs_outdir_clust)

    # Run diamond blastx on plant virus palmprint database
    blast_outdir_plantvir = os.path.join(overview_outdir,"annotate","dblastx_virpalmprints")
    if not os.path.exists(blast_outdir_plantvir):
        os.makedirs(blast_outdir_plantvir)
    blast_outpath_plantvir = os.path.join(blast_outdir_plantvir, f"{os.path.split(master_contigs)[1].split('.')[0]}_dblastx_plantvir.m8")
    blast_outpath_plantvir = run_simple_dblastx(blast_outpath_plantvir, master_contigs, PLANT_VIRUS_DB)

    # Run palmscan
    palmscan_outdir = os.path.join(overview_outdir,"annotate","palmscan")
    palmscan_outpath = run_palmscan(os.path.split(master_contigs)[1].split('.')[0], palmscan_outdir, master_contigs)

    # Run diamond blastx on palmdb
    blast_outdir_palmdb = os.path.join(overview_outdir,"annotate","dblastx_palmdb")
    if not os.path.exists(blast_outdir_palmdb):
        os.makedirs(blast_outdir_palmdb)
    blast_outpath_palmdb = os.path.join(blast_outdir_palmdb, f"{os.path.split(master_contigs)[1].split('.')[0]}_dblastx_palmdb.m8")
    blast_outpath_palmdb = run_simple_dblastx(blast_outpath_palmdb, master_contigs, PALM_DB)

    # Run diamond blastx on blast_nr
    blast_outdir_nr = os.path.join(overview_outdir,"annotate","dblastx_nr")
    if not os.path.exists(blast_outdir_nr):
        os.makedirs(blast_outdir_nr)
    blast_outpath_nr = os.path.join(blast_outdir_nr, f"{os.path.split(master_contigs)[1].split('.')[0]}_dblastx_nr.m8")
    blast_outpath_nr = run_dblasx(blast_outpath_nr, master_contigs, BLAST_NR_DB)
    filtered_blastx_fn = os.path.join(blast_outdir_nr, f"{os.path.split(master_contigs)[1].split('.')[0]}_dblastx_nr_filtered.tsv")
    filtered_blastx_fn = filter_blastx_nr(blast_outpath_nr, filtered_blastx_fn, coverage = 50, seq_id_thres= 0.85)
    # Run mmseqs search embryophyta
    mmseqs_outdir_embryo = os.path.join(overview_outdir,"annotate","mmseqs_search_embryophyta")
    if not os.path.exists(mmseqs_outdir_embryo):
        os.makedirs(mmseqs_outdir_embryo)
    mmseqs_outpath_embryo = mmseqs_easy_search_simple(master_contigs, EMBRYOPHYTA_DB, "embryophyta", mmseqs_outdir_embryo)
    # Run genomad
    genomad_outdir = os.path.join(overview_outdir,"annotate","genomad")
    genomad_outpath = run_genomad(master_contigs,GENOMAD_DB,genomad_outdir)
    # Run mmseqs search blast_nt
    mmseqs_outdir_nt = os.path.join(overview_outdir,"annotate","mmseqs_search_nt")
    if not os.path.exists(mmseqs_outdir_nt):
        os.makedirs(mmseqs_outdir_nt)
    mmseqs_outpath_nt = mmseqs_easy_search_tax(master_contigs, MMSEQS_NT, "mmseqs_nt", mmseqs_outdir_nt)

    return mmseqs_outpath_clust,blast_outpath_plantvir, palmscan_outpath, blast_outpath_palmdb, blast_outpath_nr, mmseqs_outpath_embryo, genomad_outpath, mmseqs_outpath_nt


def parse_palmdb_blastx(palmdb_blastx_fn):
    """
    Parse the dblastx output from PalmDB, store in dict
    :param palmdb_blastx_fn: path to palmdb blastx output
    :return: palmdb_blastx_dict: dict with contig as key and list of rest of blastx columns as values
    """
    palmdb_blastx_dict = {}
    with open(palmdb_blastx_fn, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] not in palmdb_blastx_dict:
                palmdb_blastx_dict[line[0]] = [line[1:]]
            else:
                palmdb_blastx_dict[line[0]].append(line[1:])

    return palmdb_blastx_dict


def filtering_hits(palmdb_blastx_dict, bitscore_thres=100, seq_id_thres= 50):
    """
    Filter hits based on bitscore and seq_id
    :param max_bitscore_hits: dict with contig as key and list of rest of blastx columns as values
    :return: filtered_dict: Dict of lists, key:contig_name, list: blastx columns of hits passing the filters
    """

    filtered_dict = {}

    for contig,hits in palmdb_blastx_dict.items():
        for hit in hits:
            if float(hit[1]) >= seq_id_thres and float(hit[10]) >= bitscore_thres:
                if contig not in filtered_dict:
                    filtered_dict[contig] = [hit]
                else:
                    filtered_dict[contig].append(hit)

    return filtered_dict

def parse_filtered_annotations(filtered_dict):
    """

    :param filtered_dict: Dict of lists,
    :return:
    """
    palmdb_filtered_dict = {}

    for contig, hits in filtered_dict.items():
        annot_list = []
        for hit in hits:
            annotation = hit[0].split("_")[0]
            annot_list.append(annotation)
        annot_set = set(annot_list)
        if len(annot_set) == 1:
            max_bitscore_hit = max(hits, key=lambda x: float(x[10]))
            hit_name = max_bitscore_hit[0]
            pid = max_bitscore_hit[1]
            bitscore = max_bitscore_hit[10]
            palmdb_string = f"{hit_name}|bs:{bitscore}|pid:{pid}"
            palmdb_filtered_dict[contig] = palmdb_string
        else:
            per_annot_dict = {s: [] for s in annot_set}
            for hit in hits:
                annotation = hit[0].split("_")[0]
                per_annot_dict[annotation].append(hit)

            annot_list = []
            for annot, hits in per_annot_dict.items():
                max_bitscore_hit = max(hits, key=lambda x: float(x[10]))
                hit_name = max_bitscore_hit[0]
                pid = max_bitscore_hit[1]
                bitscore = max_bitscore_hit[10]
                palmdb_string = f"{hit_name}|bs:{bitscore}|pid:{pid}"
                annot_list.append(palmdb_string)
            palmdb_filtered_dict[contig] = ",".join(annot_list)

    return palmdb_filtered_dict




def parse_palmscan_contigs(palm_out):
    """
    Parses the palmscan output file and returns a set of contig names.

    :param palm_out: Path to the palmscan output file with the .fev extension
    :type palm_out: str
    :return: query_set: Set of contig names extracted from the palmscan output file
    """
    palmscan_list = []

    with open(palm_out, 'r') as handle:
        for line in handle:
            query = line.split('\t')[1]
            contig = query.split("=")[1]
            palmscan_list.append(contig)

    return palmscan_list


def parse_genomad(genomad_fn):
    """
    Parses the genomad output file and returns a dict with contig as key and list [n_genes, virus_score, taxonomy] as value
    :param genomad_fn: path to genomad output file
    :return: genomad_dict: dict with contig as key and list [n_genes, virus_score, taxonomy] as value
    """
    genomad_dict = {}
    with open(genomad_fn, 'r') as f:
        for line in f:
            if line.startswith('seq_name'):
                continue

            line = line.strip().split('\t')
            n_genes = line[4]
            virus_score = line[6]
            taxonomy = line[10]
            if line[0] not in genomad_dict:
                genomad_dict[line[0]] = [n_genes, virus_score, taxonomy]
            else:
                genomad_dict[line[0]].append([n_genes, virus_score, taxonomy])

    return genomad_dict

def parse_dblastx_nr(dblasx_nr_fn):
    """
    Parses the dblastx_nr output file and returns a dict with contig as key and list of rest of blastx columns as values
    :param dblasx_nr_fn: path to dblastx_nr output file
    :return: dict with contig as key and list of rest of blastx columns as values
    """

    dblasx_nr = {}

    with open(dblasx_nr_fn, 'r') as f:
        for line in f:
            line = line.strip().split('\t')

            if line[0] not in dblasx_nr:
                dblasx_nr[line[0]] = [line[1:]]
            else:
                dblasx_nr[line[0]].append(line[1:])
    return dblasx_nr


def parse_sskingdoms(dblastx_nr_fn):
    """
    Parses the dblastx_nr output file and returns a dict with contig as key and top 5 hits sskingdoms as values
    :param dblastx_nr_fn: dict with contig as key and list of rest of blastx columns as values
    :return: dict with contig as key and top 5 hits sskingdoms as values
    """
    contig_sskingdoms = {}

    with open(dblastx_nr_fn, 'r') as f:
        for line in f:
            line = line.strip().split('\t')

            if line[0] not in contig_sskingdoms:
                contig_sskingdoms[line[0]] = [line[14]]
            else:
                contig_sskingdoms[line[0]].append(line[14])

    return contig_sskingdoms


def parse_skingdoms(dblastx_nr_fn):
    """
    Parses the dblastx_nr output file and returns a dict with contig as key and top 5 hits skingdoms as values
    :param dblastx_nr_fn: dict with contig as key and list of rest of blastx columns as values
    :return: dict with contig as key and top 5 hits skingdoms as values
    """
    contig_skingdoms = {}

    with open(dblastx_nr_fn, 'r') as f:
        for line in f:
            line = line.strip().split('\t')

            if line[0] not in contig_skingdoms:
                contig_skingdoms[line[0]] = [line[15]]
            else:
                contig_skingdoms[line[0]].append(line[15])
    return contig_skingdoms

def parse_best_hit_blastx(dblasx_nr_dict):
    """
    Parses the dblastx_nr_dict and returns a dict with contig as key and best hit as value
    :param dblasx_nr_dict:
    :return:
    """
    best_hit_blastx_dict = {}
    for contig, hits in dblasx_nr_dict.items():
        min_eval_hit = min(hits, key=lambda x: float(x[9]))
        best_hit_blastx_dict[contig] = min_eval_hit

    return best_hit_blastx_dict


def parse_embryophyta_blastn(embryophyta_blastn_fn):
    """
    Parses the embryophyta_blastn_fn and returns a dict with contig as key and list of rest of blastn columns as values
    :param embryophyta_blastn_fn:
    :return:
    """
    emb_hits_dict = {}
    with open(embryophyta_blastn_fn, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] not in emb_hits_dict:
                emb_hits_dict[line[0]] = [line[1:]]
            else:
                emb_hits_dict[line[0]].append(line[1:])

    return emb_hits_dict

def parse_pid_bitscore_embryophyta(emb_hits_dict):

    embryo_dict = {}

    for contig, hits in emb_hits_dict.items():
        min_eval_hit = min(hits, key=lambda x: float(x[9]))
        bitscore = min_eval_hit[10]
        pid = min_eval_hit[1]
        eval = min_eval_hit[9]
        query_cov = min_eval_hit[11]
        embryo_string = f"bs:{bitscore}|eval:{eval}|pid:{pid}|query_cov:{query_cov}"
        embryo_dict[contig] = embryo_string

    return embryo_dict

def parse_highscore_hmm(highest_scoring_hmm_fn):

    hmm_dict = {}
    with open(highest_scoring_hmm_fn, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            sample = line[30]
            contig = line[31]
            sample_contig = f"{sample}_{contig}"

            norm_bitscore_profile = float(line[23])
            profile_cov = float(line[28])
            contig_cov = float(line[29])
            eval = line[6]
            bitscore = line[7]
            hmm_string = f"eval:{eval}|bitscore:{bitscore}|norm_bitscore:{norm_bitscore_profile:.3f}|profile_cov:{profile_cov:.3f}|contig_cov:{contig_cov:.3f}"
            hmm_dict[sample_contig] = hmm_string
    return hmm_dict


def parse_mmseqs_search_nt(mmseqs_outpath):


    mmseqs_nt_dict = {}
    with open(mmseqs_outpath, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] not in mmseqs_nt_dict:
                mmseqs_nt_dict[line[0]] = [line[1:]]
            else:
                mmseqs_nt_dict[line[0]].append(line[1:])

    processed_mmseqs_nt_dict = {}
    for contig, hits in mmseqs_nt_dict.items():
        min_eval_hit = min(hits, key=lambda x: float(x[9]))
        bitscore = min_eval_hit[10]
        pid = min_eval_hit[1]
        qcov = min_eval_hit[11]
        lineage = min_eval_hit[15]
        eval = min_eval_hit[9]
        mmseqs_nt_string = f"bitscore:{bitscore}|pid:{pid}|qcov:{qcov}|eval:{eval}|lineage:{lineage}"
        processed_mmseqs_nt_dict[contig] = mmseqs_nt_string

    return processed_mmseqs_nt_dict





def write_annotation_file( genomad_dict, cluster_fn,palmdb_filtered_dict , palmscan_contigs, contig_skingdoms_dict,
                           contig_sskingdoms_dict,best_hit_blastx_dict, embryo_dict,palmdb_out_dict, hmm_dict, out_fn,processed_mmseqs_nt_dict):

    with open(out_fn, 'w') as out_handle, open(cluster_fn,'r') as in_handle:
        out_handle.write("Cluster\tContig\tLength\tGenomad_n_genes\tGenomad_virus_score\tGenomad_taxonomy\tPalmdb_hit(bitscore >=100, pid >=50%)\tPalmdb_taxonomy\tPlant_palmprintsDB((bitscore >=100, pid >=50%))\t"
                         "PalmScan positive\tEmbryophyta positive\tTop5_blastx_nr_sskingdoms\tTop5_blastx_nr_skingdoms\tBest_blastx_nr_acc\t"
                         "Best_blastx_nr_taxonomy\tBest_blastx_nr_taxid\tBest_blastx_nr_title\tBest_blastx_nr_bitscore|pid\thmm|eval|bs|norm_bs|prof_cov|cont_cov\tmmseqs_nt\n")
        for line in in_handle:
            line= line.strip().split('\t')
            length = line[1].split('_')[4]
            if line[1] in genomad_dict:
                n_genes = genomad_dict[line[1]][0]
                virus_score = genomad_dict[line[1]][1]
                taxonomy = genomad_dict[line[1]][2]
            else:
                n_genes = 'NA'
                virus_score = 'NA'
                taxonomy = 'NA'

            if line[1] in palmdb_filtered_dict:
                palmdb_string = palmdb_filtered_dict[line[1]]
            else:
                palmdb_string = 'NA'

            if line[1] in palmdb_out_dict:
                palmdb_out_string = palmdb_out_dict[line[1]]
                palmdb_bit_pid, palmdb_taxonomy = palmdb_out_string.rsplit("|",1)
            else:
                palmdb_bit_pid = 'NA'
                palmdb_taxonomy = 'NA'

            if line[1] in palmscan_contigs:
                palmscan_bool = 'True'
            else:
                palmscan_bool = 'False'

            if line[1] in embryo_dict:
                emb_string = embryo_dict[line[1]]
            else:
                emb_string = 'NA'

            if line[1] in contig_skingdoms_dict:
                contig_skingdoms = ",".join(contig_skingdoms_dict[line[1]])
            else:
                contig_skingdoms = 'NA'

            if line[1] in contig_sskingdoms_dict:
                contig_sskingdoms = ",".join(contig_sskingdoms_dict[line[1]])
            else:
                contig_sskingdoms = 'NA'

            if line[1] in best_hit_blastx_dict:
                best_hit_acc = best_hit_blastx_dict[line[1]][0]
                best_hit_bitscore = best_hit_blastx_dict[line[1]][10]
                best_hit_pid = best_hit_blastx_dict[line[1]][1]
                best_hit_cov = best_hit_blastx_dict[line[1]][17]
                best_hit_bs_pid = f"bs:{best_hit_bitscore}|pid:{best_hit_pid}|best_hit_cov:{best_hit_cov}"
                best_hit_sskingdom = best_hit_blastx_dict[line[1]][13]
                best_hit_skingdom = best_hit_blastx_dict[line[1]][14]
                best_hit_sphylums = best_hit_blastx_dict[line[1]][15]
                best_hit_sscinames = best_hit_blastx_dict[line[1]][12]
                best_hit_tax_id = best_hit_blastx_dict[line[1]][11]
                best_hit_taxonomy = f"{best_hit_sskingdom};{best_hit_skingdom};{best_hit_sphylums};{best_hit_sscinames}"
                best_hit_title = best_hit_blastx_dict[line[1]][16]
            else:
                best_hit_acc = 'NA'
                best_hit_bs_pid = 'NA'
                best_hit_taxonomy = 'NA'
                best_hit_title = 'NA'
                best_hit_tax_id = 'NA'


            if line[1] in hmm_dict.keys():
                hmm_string = hmm_dict[line[1]]
            else:
                hmm_string = 'NA'
            if line[1] in processed_mmseqs_nt_dict.keys():
                mmseqs_nt_string = processed_mmseqs_nt_dict[line[1]]
            else:
                mmseqs_nt_string = 'NA'

            out_handle.write(f"{line[0]}\t{line[1]}\t{length}\t{n_genes}\t{virus_score}\t{taxonomy}\t{palmdb_bit_pid}\t{palmdb_taxonomy}\t{palmdb_string}\t"
                             f"{palmscan_bool}\t{emb_string}\t{contig_sskingdoms}\t{contig_skingdoms}\t{best_hit_acc}"
                             f"\t{best_hit_taxonomy}\t{best_hit_tax_id}\t{best_hit_title}\t{best_hit_bs_pid}\t{hmm_string}\t{mmseqs_nt_string}\n")

    return out_fn

def remove_controls(annotation_fn):
    with open(annotation_fn, 'r') as f, open("genomad_palmdb_annotation_no_controls.tsv",'w') as out_handle, open("genomad_palmdb_annotation_controls.tsv",'w') as controls_handle:
        for line in f:
            if line.startswith('contig'):
                out_handle.write(line)
                controls_handle.write(line)
                continue
            new_line = line.strip().split('\t')
            if new_line[1].startswith('Neg') or new_line[1].startswith('Pos') or new_line[0].startswith('Neg') or new_line[0].startswith('Pos'):
                controls_handle.write(line)
                continue
            else:
                out_handle.write(line)

    return "genomad_palmdb_annotation_no_controls.tsv"


def parse_representatives(annotation_fn):
    representatives_set = set()

    with open(annotation_fn,'r') as f, open("representatives_annotation.txt",'w') as out_handle:
        for line in f:
            if line.startswith('contig'):
                title_line = line
                continue
            line = line.strip().split('\t')
            representatives_set.add(line[0])
        out_handle.write(title_line)
        f.seek(0)
        for line_ in f:
            if line_.startswith('contig'):
                continue
            new_line = line_.strip().split('\t')
            if new_line[1] in representatives_set:

                out_handle.write(line_)

def palmdb_parse_taxonomy(palmdb_unique_tax_fn):
    palmdb_tax_dict = {}
    with open(palmdb_unique_tax_fn, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0] == 'Label':
                continue
            if line[0] not in palmdb_tax_dict:
                palmdb_tax_dict[line[0]] = [";".join(line[1:])]
            else:
                palmdb_tax_dict[line[0]].append(";".join(line[1:]))

    return palmdb_tax_dict


def parse_pid_bitscore_tax_palmdb(palmdb_dict, palmdb_tax_dict):

    palmdb_out_dict = {}

    for contig, hits in palmdb_dict.items():
        max_bs_hit = min(hits, key=lambda x: float(x[9]))
        bitscore = max_bs_hit[10]
        pid = max_bs_hit[1]
        target_cov = max_bs_hit[13]
        palmdb_tax = max_bs_hit[0]
        taxonomy = palmdb_tax_dict[palmdb_tax]
        palmdb_out_string = f"bs:{bitscore}|pid:{pid}|tcov:{target_cov}|{taxonomy[0]}"
        palmdb_out_dict[contig] = palmdb_out_string

    return palmdb_out_dict

def filter_blastx_nr(blastx_nr_fn, blastx_nr_filtered_fn, coverage = 50, seq_id_thres= 85):
    with open(blastx_nr_fn, 'r') as f, open(blastx_nr_filtered_fn, 'w') as out_handle:
        for line in f:
            line = line.strip().split('\t')
            if float(line[18]) >= coverage and float(line[2]) >= seq_id_thres and line[14] != 'Viruses':
                out_handle.write('\t'.join(line)+'\n')

    return blastx_nr_filtered_fn

def filter_blastn_embryophyta(blastn_embryophyta_fn, blastn_embryophyta_filtered_fn, coverage = 50, seq_id_thres= 0.85):
    pass

def wrap_annotation_file(overview_outdir, mmseqs_outpath_clust,blast_outpath_plantvir, palmscan_outpath, blast_outpath_palmdb, blast_outpath_nr, mmseqs_outpath_embryo, genomad_outpath,highscore_hmm_fn,mmseqs_nt_tax_fn):
    # Process plant palmprint blastx output--check

    plantvir_blastx_dict = parse_palmdb_blastx(blast_outpath_plantvir)
    filtered_dict = filtering_hits(plantvir_blastx_dict)
    plantvir_filtered_dict = parse_filtered_annotations(filtered_dict)

    # Process palmscan output--check
    palmscan_contigs = parse_palmscan_contigs(palmscan_outpath)

    # Process genomad output
    genomad_dict = parse_genomad(genomad_outpath)

    # Process dblastx_nr output--check
    dblasx_nr_dict = parse_dblastx_nr(blast_outpath_nr)
    contig_ss_kingdoms_dict = parse_sskingdoms(blast_outpath_nr)
    contig_skingdoms_dict = parse_skingdoms(blast_outpath_nr)
    best_hit_blastx_dict = parse_best_hit_blastx(dblasx_nr_dict)

    # Process embryophyta blastn output--check
    emb_hits_dict = parse_embryophyta_blastn(mmseqs_outpath_embryo)
    embryo_dict = parse_pid_bitscore_embryophyta(emb_hits_dict)

    # Process cluster file
    cluster_fn = mmseqs_outpath_clust

    # process palmdb--check
    palmdb_blastx_dict = parse_palmdb_blastx(blast_outpath_palmdb)
    filtered_dict = filtering_hits(palmdb_blastx_dict)

    # Process palmdb unique taxonomy file
    palmdb_tax_dict = palmdb_parse_taxonomy(PALM_DB_TAX)
    palmdb_out_dict = parse_pid_bitscore_tax_palmdb(filtered_dict, palmdb_tax_dict)

    # Process hmm output
    hmm_dict = parse_highscore_hmm(highscore_hmm_fn)

    # Process mmseqs_nt output
    processed_mmseqs_nt_dict = parse_mmseqs_search_nt(mmseqs_nt_tax_fn)

    # Write annotation file
    annotated_fn = os.path.join(overview_outdir,"annotate","master_contigs.annotated.tsv")
    annotated_fn = write_annotation_file(genomad_dict, cluster_fn, plantvir_filtered_dict, palmscan_contigs, contig_skingdoms_dict,
                          contig_ss_kingdoms_dict, best_hit_blastx_dict, embryo_dict,palmdb_out_dict, hmm_dict, annotated_fn,processed_mmseqs_nt_dict)

    return annotated_fn








