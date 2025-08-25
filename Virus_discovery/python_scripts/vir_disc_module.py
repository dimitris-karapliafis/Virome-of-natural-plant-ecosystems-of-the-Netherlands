'''
Author: Dimitris Karapliafis
Date: 2023-11-01
Description:

This script is a pipeline for viral discovery from metagenomic data.
It takes as input a directory containing subdirectories with contigs.fasta files and performs the following steps:
    1.  Extracts contigs longer than 400bp using seqkit
    2.  Translates contigs to amino acid sequences using transeq
    3.  Searches for RdRp and RVMT domains using hmmsearch
    4.  Extracts contigs containing RdRp and RVMT domains
    5.  Clusters contigs using mmseqs
    6.  Extracts representative sequences from clusters
    7.  Creates a master contig file containing representative sequences

Dependencies:
    1.  seqkit
    2.  transeq
    3.  hmmsearch
    4.  mmseqs
    5.  diamond

Database dependencies:
        1.  RdRp-scan database
        2.  RVMT database
        3.  NCBI nr database
Usage:
'''

##########        IMPORTS     ##########
from sys import argv
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import upsetplot
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from rich.console import Console
import pyhmmer
import argparse
import logging
import time
from logger import logger

##########        DATABASES     ##########

hmm_rdrp = "PHD/DBs/RDRP-scan/RdRp_HMM_profile_CLUSTALO.db"
hmm_neordrp = "PHD/DBs/neo-rdrp/NeoRdRp-HMM.v1.1.hmm"
hmm_firth = "PHD/DBs/firth/conc_prof.hmm"
hmm_rvmt = "PHD/DBs/rvmt/RVMT.hmm"
diam_rdrpscan = "PHD/DBs/RDRP-scan/RdRp-scan_0.90.dmnd"
diam_neordrp = "PHD/DBs/neo-rdrp/NeoRdRp-seq.v1.1.dmnd"
genomad_db = "PHD/DBs/genomad_db"
dblast_nr = 'PHD/DBs/blast_nr/diamond_nr.dmnd'



##########        CLASSES     ##########


class Logger:
   def __init__(self, log_file):
       self.console = Console()
       self.log_file = log_file
       self.logger = logging.getLogger('Logger')
       self.logger.setLevel(logging.INFO)
       handler = logging.FileHandler(self.log_file)
       handler.setLevel(logging.INFO)
       formatter = logging.Formatter('%(asctime)s - %(message)s')
       handler.setFormatter(formatter)
       self.logger.addHandler(handler)

   def log(self, message):
       self.console.log(message)
       self.logger.info(message)

   def start_timer(self):
       self.start_time = time.time()

   def stop_timer(self):
       end_time = time.time()
       execution_time = end_time - self.start_time
       self.log(f"Execution time: {execution_time} seconds")

class hmmsearch_parser:
    """
    Class for parsing hmmsearch output files.

    Attributes:
        data (dict): A dictionary containing the parsed data from the hmmscan output file.
        hmm_output_file (str): Path to the hmmscan output file.

    Methods:
        parse_output(hmm_output_file): Parses the hmmsearch output file and returns a dictionary.
        calculate_coverage(data): Calculates the coverage of all domains in a profile.
        get_contig(contig_name): Returns all profiles and domains for a given contig.
        export_processed_file(data, outfile, p_cov_threshold=0): Exports the processed hmmscan output file.
    """

    def __init__(self, hmm_raw, hmm_processed):
        """
        Constructor for the hmmsearch_parser class.

        :param hmm_raw: Path to the raw hmmsearch output file.
        :type hmm_raw: str
        :param hmm_processed: Path to the processed output file.
        :type hmm_processed: str
        """
        self.data = {}
        self.hmm_output_file = hmm_raw
        parsed_data = self.parse_output(self.hmm_output_file)
        parsed_data = self.calculate_norm_bitscore_profile(parsed_data)
        parsed_data = self.calculate_norm_bitscore_contig(parsed_data)
        parsed_data = self.calculate_norm_bitscore_custom(parsed_data)
        self.data = self.calculate_coverage(parsed_data)

        self.export_processed_file(self.data, hmm_processed)

    def parse_output(self, hmm_raw_out):
        """
        Parse hmmsearch output file and return a dictionary with the following structure:
        {contig_name: {profile_name: [[list of domain data]]}}

       :param hmm_output_file: Path to the hmmscan output file.
       :type hmm_output_file: str
       :return: Dictionary with parsed data.
       :rtype: dict
       """

        with open(hmm_raw_out, 'r') as hmm_file:

            for line in hmm_file:
                if line.startswith('#') or line.startswith('t_name'):
                    continue

                tmp_line = line.strip().split()
                for element in tmp_line:
                    element.replace(' ', '')

                desc = tmp_line[22:]
                tmp_line = tmp_line[:22] + [' '.join(desc)]

                if tmp_line[0] not in self.data:
                    self.data[tmp_line[0]] = {tmp_line[3]: [[col for col in tmp_line[1:]]]}
                else:
                    if tmp_line[3] not in self.data[tmp_line[0]]:
                        self.data[tmp_line[0]][tmp_line[3]] = [[col for col in tmp_line[1:]]]
                    else:
                        self.data[tmp_line[0]][tmp_line[3]].append([col for col in tmp_line[1:]])

        return self.data

    def calculate_norm_bitscore_profile(self, data):
        """
        Calculates the BitScore/Length for all domains in a profile.
        :param data:
        :return:
        """

        for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    for domain in domains:
                        model_length = float(domain[4])
                        domain_bitscore = float(domain[6])
                        norm_bitscore = domain_bitscore/model_length
                        domain.append(norm_bitscore)

        return data

    def calculate_norm_bitscore_contig(self, data):
        """
        Calculates the BitScore/Length for the contig size.
        :param data:
        :return:
        """

        for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    for domain in domains:
                        contig_length = float(domain[1])
                        domain_bitscore = float(domain[6])
                        norm_bitscore = domain_bitscore/contig_length
                        domain.append(norm_bitscore)

        return data

    def calculate_norm_bitscore_custom(self, data):
        """
        Calculates the 2*BitScore/Length of contig + Length of profile for all domains in a profile.
        :param data:
        :return:
        """

        for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    for domain in domains:
                        model_length = float(domain[4])
                        contig_length = float(domain[1])
                        domain_bitscore = float(domain[6])
                        norm_bitscore = (domain_bitscore/model_length) + (domain_bitscore/contig_length)
                        domain.append(norm_bitscore)

        return data

    def calculate_coverage(self, data):
        """
        Calculates the % coverage of all domains in a profile.

        :param data: Dictionary with parsed data.
        :type data: dict
        :return: Dictionary with parsed data and domain coverage.
        :rtype: dict
        """
        # TODO: This needs to be thoroughly tested, the %perc coverage is known to be difficult to calculate due to
        # TODO: many diffent scenarios of overlap between domains. See this issue for more info:
        # TODO: COntig coverage takes into account the hmm from-to positions, not the ali_from-to positions, so it can be
        # TODO: greater than 1. See this issue for more info:
        # https://github.com/althonos/pyhmmer/issues/27
        overlap_set = set()
        cov_set = set()
        for contig, profiles in data.items():

            for profile, domains in profiles.items():


                query_length = int(domains[0][4])
                contig_length = int(domains[0][1])
                bitscore = float(domains[0][6])
                for i in range(len(domains)):  # iterate over domains
                    prev_dom_start = int(domains[i][14])
                    prev_dom_end = int(domains[i][15])
                    cov_range = list(range(prev_dom_start, prev_dom_end + 1))
                    for pos in cov_range:
                        cov_set.add(pos)

                    if i == len(domains)-1:
                        domain_coverage = len(cov_set) / query_length
                        contig_coverage = len(cov_set) / contig_length
                        domains[i].append(bitscore/len(cov_set))
                        domains[i].append(len(cov_set))
                        domains[i].append(domain_coverage)
                        domains[i].append(contig_coverage)

                    cov_set = set()

        return data  # return data with domain coverage added

    def get_contig(self, contig_name):
        """
        Returns all profiles and domains for a given contig.

        :param contig_name: Name of the contig.
        :type contig_name: str
        :return: Dictionary with profiles and domains.
        :rtype: dict
        """
        return self.data.get(contig_name, {})

    def export_processed_file(self, data, outfile, p_cov_threshold=0):
        """
        Exports the processed hmmsearch output file.

        :param data: Dictionary with parsed data.
        :type data: dict
        :param outfile: Path to the output file.
        :type outfile: str
        :param p_cov_threshold: Minimum profile coverage threshold, defaults to 0.
        :type p_cov_threshold: int, optional
        :return: Path to the output file.
        :rtype: str
        """
        title_line = ["t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                      "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                      "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                      "description of target", "norm_bitscore_profile", "norm_bitscore_contig", 'norm_bitscore_custom',
                      "ID_score",'aln_length','profile_coverage', "contig_coverage"]

        line_list = []
        with open(outfile, 'w') as out:
            for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    if len(profiles[profile]) > 1:
                        line_list.append([contig] + domains[-1])
                    else:
                        for domain in domains:
                            line_list.append([contig] + domain)
            out.write("\t".join(title_line) + '\n')

            for line in line_list:
                if line[-1] > p_cov_threshold:
                    j_line = "\t".join(str(l) for l in line)
                    out.write(j_line + '\n')
        return outfile


##########        FUNCTIONS     ##########



def extract_dir_names(in_path):
    """
    Extracts directory names from a given path

    :param in_path:
    :type in_path: str
    :return: list of directory names [str]
    :rtype: list
    """
    dir_names = []

    for content in os.listdir(in_path):
        dir_names.append(content)

    return dir_names


def create_outdir(in_path, dir_names):
    """
    Creates output directories

    :param in_path:
    :type in_path: str
    :param dir_names:
    :type dir_names: list
    :return: list of output paths
    :rtype: list

    """
    parent_dir = os.path.split(in_path)[0]
    out_paths = []

    for filename in dir_names:

        out_path = os.path.join(parent_dir, f"results/{filename}_results")

        if not os.path.exists(out_path):
            os.makedirs(out_path)
        out_paths.append(out_path)

    return out_paths


def run_seqkit_filt(seqkit_out_path, contig_path, dir_name):
    """
    Runs seqkit to filter contigs shorter than 400bp

    :param seqkit_out_path: Path to the filtered contig file
    :type seqkit_out_path: str
    :param contig_path: Path to the contig file
    :type contig_path: str
    :param dir_name: Name of the directory-Sample
    :type dir_name: str
    :return: seqkit_out_path: Path to the filtered contig file
    :rtype: str
    """

    if not os.path.exists(seqkit_out_path):
        try:
            cmd_seqkit = f"seqkit seq -m 400 {contig_path} -o {seqkit_out_path}"
            res_map = subprocess.check_output(cmd_seqkit,
                                              shell=True)
        except(subprocess.CalledProcessError):
           logger.log(
                f"Something went wrong with seqkit for sample:{os.path.split(dir_name)[0]}")

    return seqkit_out_path


def run_transeq(transeq_out_path, contig_path, dir_name, gen_code=1):
    """
    Runs transeq to translate contigs to amino acid sequences in all 6 frames for a given genetic code

    :param transeq_outpath: Path to the translated contig file
    :type transeq_outpath: str
    :param contig_path: Path to the contig file
    :type contig_path: str
    :param dir_name: Name of the directory-Sample
    :type dir_name: str
    :param gen_code: A supported genetic code for transeq, defaults to 1
    :type gen_code: int, optional
    :return: transeq_out_path: Path to the translated contig file
    """

    if not os.path.exists(transeq_out_path):
        try:
            cmd_transeq = f"transeq {contig_path} {transeq_out_path} -table {gen_code} -frame 6 -clean"
            res_map = subprocess.check_output(cmd_transeq,
                                              shell=True)
        except(subprocess.CalledProcessError):
           logger.log(
                f"Something went wrong with transeq for sample:{os.path.split(dir_name)[0]}")

    return transeq_out_path


def run_pyhmmer_search(hmmsearch_out_path, seq_file, hmm_file):
    """
    Runs pyhmmer hmmsearch on a given sequence file and hmm profile database

    :param hmmsearch_out: Path to store the hmmsearch output
    :param seq_file: Path to the sequence file
    :param hmm_file: Path to the hmm profile database
    :return: hmmsearch_out: Path to the hmmsearch output

    TODO: 1. Add option to run hmmsearch on long sequences (longer than 100kb) as pyhmmer.Pipeline is not able to handle
    TODO: long sequences.  See: https://pyhmmer.readthedocs.io/en/latest/api/plan7.html#pyhmmer.plan7.LongTargetsPipeline
    TODO: 2. Parameters are now hardcoded, add option to change them
    """

    if not os.path.exists(hmmsearch_out_path):

        with pyhmmer.plan7.HMMPressedFile(hmm_file) as handle:
            hmms = list(handle)

        with pyhmmer.easel.SequenceFile(seq_file, digital=True) as handle:
            db = list(handle)

        with open(hmmsearch_out_path, 'wb') as handle:
            title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                          "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                          "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                          "description of target"]
            handle.write("\t".join(title_line).encode("utf-8") + b"\n")

            for result in pyhmmer.hmmer.hmmsearch(hmms, db, cpus=20, E=1e-5, incdomE=1e-5, domE=1e-5, incE=1e-5,
                                                  Z=1000000):
                result.write(handle, format="domains", header=False)

    return hmmsearch_out_path


def wrap_hmmsearch(transeq_out_path, rvmt_db, rdrpscan_db, hmm_outdir, sample, analysis_dict):
    """
    Wrapper function for running hmmsearch on a given sample for both RVMT and RdRp-Scan pHMM databases.

    :param transeq_out_path: Path to the translated contig file
    :type transeq_out_path: str
    :param rvmt_db: Path to the RVMT database
    :type rvmt_db: str
    :param rdrpscan_db: Path to the RdRp-Scan database
    :type rdrpscan_db: str
    :param hmm_outdir: Path to the hmmsearch output directory
    :type hmm_outdir: str
    :param sample: Name of the sample
    :type sample: str
    :param analysis_dict: Dictionary containing the analysis results
    :type analysis_dict: dict
    :return: analysis_dict: Dictionary containing the analysis results
    :rtype: dict

    Analysis dict structure:
    {sample: {rvmt: {contig: {profile: [[list of domain data]]}}, rdrp: {contig: {profile: [[list of domain data]]}}}}
    This choice is made to store the full analysis data in a structure that updates with every iteration
    """

    analysis_dict[sample] = {}

    out_hmm_rvmt = os.path.join(hmm_outdir, f"{sample}_rvmt_hmmsearch.txt")

    if not os.path.exists(out_hmm_rvmt):
        run_pyhmmer_search(out_hmm_rvmt, transeq_out_path, rvmt_db)

    rvmt_out = parse_raw_hmm(out_hmm_rvmt, 'rvmt', sample)

    analysis_dict[sample]['rvmt'] = rvmt_out

    out_hmm_rdrp = os.path.join(hmm_outdir, f"{sample}_rdrp_hmmsearch.txt")

    if not os.path.exists(out_hmm_rdrp):
        run_pyhmmer_search(out_hmm_rdrp, transeq_out_path, rdrpscan_db)

    rdrp_out = parse_raw_hmm(out_hmm_rdrp, 'rdrp', sample)
    analysis_dict[sample]['rdrp'] = rdrp_out

    return analysis_dict


def parse_raw_hmm(hmm_fn_path, db_name, sample):
    """
    Parses the raw hmmsearch output file and returns a processed file, containing additional information, desc below.

    :param hmm_fn: Path to the hmmsearch output file
    :type hmm_fn: str
    :param db_name: Name of the database used for hmmsearch
    :type db_name: str
    :param sample: Name of the sample
    :type sample: str
    :return: final_hmm_path Path to the final processed hmmsearch output file
    :rtype: str

    First the class hmmsearch_parser is called to parse the raw hmmsearch output file, and add the domains coverage
    percentage as a new column. A processed file is written for inspection, as we need to optimise the class in the
    future. Additional columns are added to the processed file, containing information about the sample, the contig
    name, the database name, the contig length and the kmer coverage. The final  processed file is then used for
    downstream analysis.
    """

    hmm_processed_fn = f"{os.path.split(hmm_fn_path)[1].split('.')[0]}_processed.txt"
    hmm_processed_path = os.path.join(os.path.split(hmm_fn_path)[0], hmm_processed_fn)
    # Intermediate file, can be removed in the future

    hmmsearch_parser(hmm_fn_path, hmm_processed_path)

    final_out_fn = f"{os.path.split(hmm_fn_path)[1].split('.')[0]}_final.txt"
    final_hmm_path = os.path.join(os.path.split(hmm_fn_path)[0], final_out_fn)
    # Final file, used for downstream analysis
    with open(hmm_processed_path, 'r') as in_handle, open(final_hmm_path, 'w') as out_handle:
        title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                      "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                      "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                      "description of target",  "norm_bitscore_profile", "norm_bitscore_contig", 'norm_bitscore_custom',
                      "ID_score",'aln_length','profile_coverage', "contig_coverage", 'sample', 'cont_name', 'db_name',
                      'cont_length', 'kmer_coverage']
        out_handle.write("\t".join(title_line) + "\n")

        for line in in_handle:
            if line.startswith("#") or line.startswith("t_name"):
                continue
            else:
                no_whitespace_line = line.strip().split("\t")
                hit_name = no_whitespace_line[0]
                tmp_cont = hit_name.split("_")
                if tmp_cont[0] == "NODE":
                    cont_name = ""
                    for i in range(0, 6):
                        cont_name += tmp_cont[i] + "_"

                    cont_name = cont_name[:-1]
                    cont_length = cont_name.split("_")[3]
                    kmer_cov = cont_name.split("_")[5]
                else:
                    cont_name = hit_name
                    cont_length = "NA"
                    kmer_cov = "NA"

                no_whitespace_line.append(sample)
                no_whitespace_line.append(cont_name)
                no_whitespace_line.append(db_name)
                no_whitespace_line.append(cont_length)
                no_whitespace_line.append(kmer_cov)

                out_handle.write("\t".join(no_whitespace_line) + "\n")
                no_whitespace_line = []

    return final_hmm_path


def parse_hmm_to_set(hmm_processed_fn_path):
    """
    Extracts the contig names from the hmmsearch output file, combines to sample_contig format and returns a set.

    :param hmm_processed_fn_path: Path to the hmmsearch output file
    :type hmm_processed_fn_path: str
    :return: contigs_set: Set of unique contig names
    :rtype: set

    Example: sample: C-AM-A20, Contig: NODE_1_length_1000_cov_1.000000 -> C-AM-A20_NODE_1_length_1000_cov_1.000000
    """
    contigs_set = set()

    with open(hmm_processed_fn_path, 'r') as in_handle:
        for line in in_handle:
            if line.startswith('#t_name'):
                continue
            sample = line.split('\t')[24]
            contig = line.split('\t')[25]
            sample_contig = f"{sample}_{contig}"
            contigs_set.add(sample_contig)

    return contigs_set


def parse_palm_to_set(palm_out):
    """
    Parses the palmscan output file and returns a set of contig names.

    :param palm_out: Path to the palmscan output file with the .fev extension
    :type palm_out: str
    :return: query_set: Set of contig names extracted from the palmscan output file
    """
    query_set = set()
    with open(palm_out, 'r') as handle:
        for line in handle:
            query = line.split('\t')[1]
            contig = query.split("=")[1]
            query_set.add(contig)
    return query_set


def extract_contigs(contig_names, contig_file, sample, record_list):
    """
    Extracts contigs from a contig file based on a set of contig names and writes them to a new file.

    :param contig_names: list of contig names
    :type contig_names: list
    :param contig_file: Path to the contig file
    :type contig_file: str
    :param out_file: Path to the output file
    :type out_file: str
    :param sample: Name of the sample
    :type sample: str
    :return: out_file: Path to the output file
    """

    records = SeqIO.parse(contig_file, "fasta")

    for record in records:
        if record.id in contig_names:
            record.id = f"{sample}_{record.id}"
            record.description = ''
            record_list.append(record)

    return record_list


# def wrap_extract_contigs(in_path, general_outdir, sum_fn):
#     """
#     Wrapper function for extracting contigs from a given set of samples.
#
#     :param in_path: Path to the input directory containing the samples
#     :type in_path: str
#     :param general_outdir: Path to the general output directory
#     :type general_outdir: str
#     :param sum_fn: Path to the summary file containing the hmmsearch results
#     :type sum_fn: str
#     :return: contigs_outpath: Path to the master contig file
#     :rtype: str
#     """
#     tup_set = set()
#     cont_list = []
#     samples = extract_dir_names(in_path)
#     logger.log(f"Extracting contigs from samples: {' '.join(samples)}")
#     contigs_outpath = os.path.join(general_outdir, "master_contigs.fasta")
#     print(f"contigs_outpath: {contigs_outpath}")
#     if not os.path.exists(contigs_outpath):
#         with open(sum_fn, 'r') as handle:
#             for line in handle:
#                 if line.startswith('#'):
#                     continue
#                 line = line.strip().split('\t')
#                 sample = line[24]
#                 contig = line[25]
#                 track_tup = (contig, sample)
#                 tup_set.add(track_tup)
#         with open(contigs_outpath, 'w') as out_handle:
#             for sample in samples:
#                 in_dir = os.path.join(in_path, sample)
#                 contig_path = os.path.join(in_dir, "contigs.fasta")
#                 print(f"contig_inpath: {contig_path}")
#                 for contig, sample_name in tup_set:
#                     if sample_name == sample:
#                         cont_list.append(contig)
#
#                 extract_contigs(cont_list, contig_path, out_handle, sample)
#                 cont_list = []
#
#     return contigs_outpath


def wrap_extract_contigs(in_path, general_outdir, sum_fn):
    """
    Wrapper function for extracting contigs from a given set of samples.

    :param in_path: Path to the input directory containing the samples
    :type in_path: str
    :param general_outdir: Path to the general output directory
    :type general_outdir: str
    :param sum_fn: Path to the summary file containing the hmmsearch results
    :type sum_fn: str
    :return: contigs_outpath: Path to the master contig file
    :rtype: str
    """
    cont_dict = {}
    record_list = []


    samples = extract_dir_names(in_path)
    logger.log(f"Extracting contigs from samples: {' '.join(samples)}")
    contigs_outpath = os.path.join(general_outdir, "master_contigs.fasta")
    if not os.path.exists(contigs_outpath):
        with open(sum_fn, 'r') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                sample = line[30]
                contig = line[31]
                if sample not in cont_dict:
                    cont_dict[sample] = [contig]
                else:
                    cont_dict[sample].append(contig)

        with open(contigs_outpath, 'w') as out_handle:
            for sample in samples:
                in_dir = os.path.join(in_path, sample)
                contig_path = os.path.join(in_dir, "contigs.fasta")
                for dict_sample, cont_list in cont_dict.items():
                    if dict_sample == sample:
                        extract_contigs(cont_list, contig_path,sample, record_list)
            SeqIO.write(record_list, out_handle, "fasta")


    return contigs_outpath


def wrap_extract_trans(in_path, overview_outdir, general_outdir, sum_fn):
    """
    Wrapper function for extracting translated contigs from a given set of samples.

    :param in_path: Path to the input directory containing the samples
    :type in_path: str
    :param overview_outdir: Path to the overview output directory
    :type overview_outdir: str
    :param general_outdir: Path to the general output directory
    :type general_outdir: str
    :param sum_fn: Path to the summary file containing the hmmsearch results
    :type sum_fn: str
    :return: contigs_outpath: Path to the master contig file
    :rtype: str
    """
    tup_set = set()
    cont_list = []
    record_list = []
    samples = extract_dir_names(in_path)
    logger.log(f"Extracting translated_contigs from samples: {' '.join(samples)}")

    contigs_outpath = os.path.join(overview_outdir, "master_contigs_trans.fasta")

    if not os.path.exists(contigs_outpath):
        with open(sum_fn, 'r') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                sample = line[30]
                trans_contig = line[0]
                track_tup = (trans_contig, sample)
                tup_set.add(track_tup)
        with open(contigs_outpath, 'w') as out_handle:
            for sample in samples:
                outfile_name = f"{sample}_results"
                trans_contig_path = os.path.join(general_outdir, outfile_name, "transeq", f"{sample}_transeq.fasta")
                for contig, sample_name in tup_set:
                    if sample_name == sample:
                        cont_list.append(contig)

                record_list =extract_contigs(cont_list, trans_contig_path, sample, record_list)
                cont_list = []
            SeqIO.write(record_list, out_handle, "fasta")

    return contigs_outpath


def parse_blastx_vir(blastx_dir_path, samples, general_outdir):
    """
    Parses the blastx virus output files and returns a set of contig names.

    :param blastx_dir_path: Path to the blastx output directory
    :type blastx_dir_path: str
    :param samples: List of sample names
    :type samples: list
    :param general_outdir: Path to the general output directory
    :type general_outdir: str

    :return: blastx_set: Set of contig names
    """

    blastx_set = set()

    for file in os.listdir(blastx_dir_path):
        if file.split("_")[0] in samples:
            blastx_path = os.path.join(blastx_dir_path, file)

            with open(blastx_path, "r") as f:
                for line in f:
                    new_line = line.strip().split("\t")
                    blastx_set.add(f"{file.split('_')[0]}_{new_line[0]}")

    return blastx_set


def extract_contigs_blastx(blastx_contig_set, blastx_contig_out, in_path):
    """
    Extracts contigs from a contig file based on a set of contig names and writes them to a new file.

    :param blastx_contig_set: Set of contig names
    :type blastx_contig_set: set
    :param blastx_contig_out: Path to the output file
    :type blastx_contig_out: str
    :param in_path: Path to the input directory containing the samples
    :type in_path: str
    :return: blastx_contig_out: Path to the output file
    :rtype: str
    """

    samp_cont_dict = {}

    for element in list(blastx_contig_set):
        tmp_line = element.strip().split('_', 1)
        sample = tmp_line[0]
        contig = tmp_line[1]

        if sample not in samp_cont_dict:
            samp_cont_dict[sample] = [contig]
        else:
            samp_cont_dict[sample].append(contig)

    # This  way, we can extract all contigs from a sample at once
    with open(blastx_contig_out, 'w') as out_handle:

        for sample, contigs in samp_cont_dict.items():
            in_dir = os.path.join(in_path, sample)
            contig_path = os.path.join(in_dir, "contigs.fasta")
            extract_contigs(contigs, contig_path, out_handle, sample)

    return blastx_contig_out


def write_summary_file(analysis_dict, outdir):
    """
    Writes a summary file from the analysis dictionary.

    :param analysis_dict: Dictionary containing the analysis results
    :type analysis_dict: dict
    :param outdir: Path to the output directory
    :type outdir: str
    :return: summ_fn: Path to the summary file
    :rtype: str
    """

    summ_fn = os.path.join(outdir, "summary.tsv")
    with open(summ_fn, 'w') as out_handle:
        title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                      "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                      "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                      "description of target", "norm_bitscore_profile", "norm_bitscore_contig", 'norm_bitscore_custom',"ID_score",'profile_coverage', "contig_coverage", 'sample', 'cont_name', 'db_name',
                      'cont_length', 'kmer_coverage']

        out_handle.write("\t".join(title_line) + "\n")

        for sample, hmm_dict in analysis_dict.items():
            for hmm_db_name, hmm_outfile in hmm_dict.items():
                with open(hmm_outfile, 'r') as handle:
                    for line in handle:
                        if line.startswith('#'):
                            continue
                        out_handle.write(line)

    return summ_fn


def write_filt_summary_file(summary_file_path, eval_thresh=1e-5):
    """
    Writes a filtered summary file based on the e-value threshold.

    :param summary_file_path: Path to the summary file
    :type summary_file_path: str
    :param eval_thresh: E-value threshold, defaults to 1e-10
    :type eval_thresh: float, optional
    :return: filter_out: Path to the filtered summary file
    """

    filter_out = os.path.join(os.path.split(summary_file_path)[0], f"summary_filtered_{eval_thresh}.tsv")

    with open(summary_file_path, 'r') as handle, open(filter_out, 'w') as out_handle:
        title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                      "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                      "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                      "description of target",  "norm_bitscore_profile", "norm_bitscore_contig", 'norm_bitscore_custom',
                      "ID_score",'aln_length','profile_coverage', "contig_coverage", 'sample', 'cont_name', 'db_name',
                      'cont_length', 'kmer_coverage']

        out_handle.write("\t".join(title_line) + "\n")

        for line in handle:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')

            if float(line[6]) <= eval_thresh:
                out_handle.write('\t'.join(line) + '\n')

    return filter_out


def highest_scoring_hmms(hmm_outfile, overview_outdir):
    """
    Writes a file containing the highest scoring hmm hits for each contig.

    :param hmm_outfile: Path to the hmm output
    :type hmm_outfile: str
    :param overview_outdir: Path to the overview output directory
    :type overview_outdir: str
    :return: out_handle_all: Path to the output file
    """

    hmm_dict = {}
    highest_scoring_list = {}
    with open(hmm_outfile) as in_handle:
        for line in in_handle:
            if line.startswith('#'):
                title_line = line
                continue
            line = line.strip().split('\t')
            trans_contig_name = line[0]
            if trans_contig_name not in hmm_dict:
                hmm_dict[trans_contig_name] = [line]
            else:
                hmm_dict[trans_contig_name].append(line)

    outfile_all = os.path.join(overview_outdir, f"master_contigs_hmm_only_top_hits.tsv")
    with open(outfile_all, 'w') as out_handle_all:
        out_handle_all.write(title_line)
        for trans_contig_name, hit_lines in hmm_dict.items():
            # Highest profile coverage selection
            max_list = max(hit_lines, key=lambda x: float(x[28]))
            out_handle_all.write('\t'.join(max_list) + '\n')

    return outfile_all


# def blastx_annotate_mmseqs(mmseqs_summary, blastx_dir_path, overview_outdir):
#
#     blastx_contig_list = []
#     with open(mmseqs_summary, 'r') as handle, open(os.path.join(overview_outdir,
#     f"master_contigs_mmseqs_annotated.tsv"), 'w') as out_handle:
#         title_line = ["#cluster_name", "cluster_member","sseqid", "pident", "length", "mismatch", "gapopen",
#         "qstart", "qend", "sstart", "send", "eval", "bitscore", "staxids", "sscinames", 'sskingdoms', 'skingdoms',
#         'sphylums', "stitle", "qcovhsp", "slen", "qlen"]
#         out_handle.write('\t'.join(title_line) + '\n')
#
#         for mmseqs_line in handle:
#             if mmseqs_line.startswith('#'):
#                 continue
#
#             line = mmseqs_line.strip().split('\t')
#             contig_line = line[1].split('_')
#
#             if len(contig_line) == 7:
#                 sample, contig = line[1].split('_', 1)
#             else:
#                 continue
#             blastx_path = os.path.join(blastx_dir_path, f"{sample}_vir.txt")
#
#             with open(blastx_path, 'r') as f:
#                 for blast_line in f:
#                     new_line = blast_line.strip().split('\t')
#                     blastx_contig = new_line[0]
#                     blastx_contig_list.append(blastx_contig)
#                     if blastx_contig == contig:
#                         final_line = mmseqs_line.strip().split('\t') + new_line[1:]
#                         out_handle.write('\t'.join(final_line) + '\n')
#                         break
#                 if contig not in blastx_contig_list:
#                     final_line = mmseqs_line.strip().split('\t') + ['NA' for i in range(0, 20)]
#                     out_handle.write('\t'.join(final_line) + '\n')
#             blastx_contig_list = []


def blastx_annotate_mmseqs(mmseqs_summary, master_blastx_path, annotate_dir):
    """
    Annotates the mmseqs output file with the blastx output file.

    :param mmseqs_summary: Path to the mmseqs output file
    :type mmseqs_summary: str
    :param master_blastx_path: Path to the blastx output file
    :type master_blastx_path: strq
    :param annotate_dir: Path to the annotation output directory
    :type annotate_dir: str
    :return: Path to the annotated mmseqs output file
    :rtype: str

    """

    blastx_contig_list = []
    with open(mmseqs_summary, 'r') as handle, open(os.path.join(annotate_dir, f"master_contigs_mmseqs_annotated.tsv"),
                                                   'w') as out_handle:

        title_line = ["#cluster_name", "cluster_member", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                      "qend", "sstart", "send", "eval", "bitscore", "staxids", "sscinames", 'sskingdoms', 'skingdoms',
                      'sphylums', "stitle", "qcovhsp", "slen", "qlen"]

        out_handle.write('\t'.join(title_line) + '\n')

        for mmseqs_line in handle:
            if mmseqs_line.startswith('#'):
                continue
            master_contig = mmseqs_line.strip().split('\t')[1]

            with open(master_blastx_path, 'r') as f:
                for blast_line in f:
                    new_line = blast_line.strip().split('\t')
                    blastx_contig = new_line[0]
                    blastx_contig_list.append(blastx_contig)
                    if blastx_contig == master_contig:
                        final_line = mmseqs_line.strip().split('\t') + new_line[1:]
                        out_handle.write('\t'.join(final_line) + '\n')
                        break

                if master_contig not in blastx_contig_list:
                    final_line = mmseqs_line.strip().split('\t') + ['NA' for i in range(0, 20)]
                    out_handle.write('\t'.join(final_line) + '\n')
            blastx_contig_list = []

    return os.path.join(annotate_dir, f"master_contigs_mmseqs_annotated.tsv")


def filtered_annotations(annotated_clust_file, annotate_dir):
    """
    Filters the annotated mmseqs output file based on the quality of the annotation: perc identity >= 80, 50 <= perc
    identity < 80, perc identity < 50.

    :param annotated_clust_file: Path to the annotated mmseqs output file
    :type annotated_clust_file: str
    :param annotate_dir: Path to the annotation output directory
    :type annotate_dir: str
    :return: None
    """

    cluster_dict = {}
    with open(annotated_clust_file, 'r') as in_handle:
        for line in in_handle:
            if line.startswith('#'):
                title_line = line
                continue
            splitted_line = line.strip().split('\t')
            cluster_name = splitted_line[0]
            cluster_annot = splitted_line[1:]
            if cluster_name not in cluster_dict:
                cluster_dict[cluster_name] = [cluster_annot]
            else:
                cluster_dict[cluster_name].append(cluster_annot)
    non_vir_cluster_dict = {}
    vir_cluster_dict = {}
    for cluster, annot_list in cluster_dict.items():
        for annot in annot_list:
            if annot[14] in ['Bacteria', 'Archaea', 'Eukaryota'] or annot[1].startswith("Pos") or annot[1].startswith(
                    "Neg"):
                non_vir_cluster_dict[cluster] = annot_list

    for cluster, annot_list in cluster_dict.items():
        if cluster not in non_vir_cluster_dict:
            vir_cluster_dict[cluster] = annot_list

    with open(os.path.join(annotate_dir, "master_contigs_mmseqs_annotated_1st_round.tsv"), 'w') as out_handle:
        out_handle.write(title_line)
        for cluster, annot_list in vir_cluster_dict.items():
            for annot in annot_list:
                out_handle.write(cluster + "\t" + '\t'.join(annot) + '\n')

    high_qual_cluster_dict = {}
    mid_qual_cluster_dict = {}
    low_qual_cluster_dict = {}
    na_qual_cluster_dict = {}
    for cluster, annot_list in vir_cluster_dict.items():
        for annot in annot_list:
            if annot[2] != 'NA':
                if float(annot[2]) >= 80:
                    high_qual_cluster_dict[cluster] = annot_list
                    break
                elif float(annot[2]) >= 50 and float(annot[2]) < 80:
                    mid_qual_cluster_dict[cluster] = annot_list
                    break
                else:
                    low_qual_cluster_dict[cluster] = annot_list
            else:
                na_qual_cluster_dict[cluster] = annot_list

    with open(os.path.join(annotate_dir, "master_contigs_mmseqs_annotated_non_vir.tsv"), 'w') as out_handle:
        out_handle.write(title_line)
        for cluster, annot_list in non_vir_cluster_dict.items():
            for annot in annot_list:
                out_handle.write(cluster + "\t" + '\t'.join(annot) + '\n')

    cluster_counter = 0
    with open(os.path.join(annotate_dir, "master_contigs_mmseqs_annotated_high_qual.tsv"), 'w') as out_handle:
        out_handle.write(title_line)
        for cluster, annot_list in high_qual_cluster_dict.items():
            cluster_counter += 1
            for annot in annot_list:
                out_handle.write(cluster + "\t" + '\t'.join(annot) + '\n')
    high_qual_no = cluster_counter
    cluster_counter = 0
    with open(os.path.join(annotate_dir, "master_contigs_mmseqs_annotated_mid_qual.tsv"), 'w') as out_handle:
        out_handle.write(title_line)

        for cluster, annot_list in mid_qual_cluster_dict.items():
            if cluster not in high_qual_cluster_dict:
                cluster_counter += 1
                for annot in annot_list:
                    out_handle.write(cluster + "\t" + '\t'.join(annot) + '\n')
    mid_qual_no = cluster_counter
    cluster_counter = 0

    with open(os.path.join(annotate_dir, "master_contigs_mmseqs_annotated_low_qual.tsv"), 'w') as out_handle:
        out_handle.write(title_line)
        for cluster, annot_list in low_qual_cluster_dict.items():
            if cluster not in high_qual_cluster_dict and cluster not in mid_qual_cluster_dict:
                cluster_counter += 1
                for annot in annot_list:
                    out_handle.write(cluster + "\t" + '\t'.join(annot) + '\n')
    low_qual_no = cluster_counter
    cluster_counter = 0

    with open(os.path.join(annotate_dir, "master_contigs_mmseqs_annotated_na_qual.tsv"), 'w') as out_handle:
        out_handle.write(title_line)
        for cluster, annot_list in na_qual_cluster_dict.items():
            if cluster not in high_qual_cluster_dict and cluster not in mid_qual_cluster_dict and cluster not in\
                    low_qual_cluster_dict:
                cluster_counter += 1
                for annot in annot_list:
                    out_handle.write(cluster + "\t" + '\t'.join(annot) + '\n')
    na_qual_no = cluster_counter
    cluster_counter = 0
    logger.log(f"Number of viral clusters: {len(vir_cluster_dict)}")
    logger.log(f"Number of non-viral clusters: {len(non_vir_cluster_dict)}")
    logger.log(f"Number of high quality annotated clusters: {high_qual_no}")
    logger.log(f"Number of mid quality annotated clusters: {mid_qual_no}")
    logger.log(f"Number of low quality annotated clusters: {low_qual_no}")
    logger.log(f"Number of NA quality annotated clusters: {na_qual_no}")


def postprocessing(master_contig_path, overview_outdir, master_blastx_path):
    """
    Postprocessing of the analysis results. Runs mmseqs clustering, annotates the clusters with blastx and filters the
     annotations based on the quality of the annotation.

    :param master_contig_path: Path to the master contig file
    :type master_contig_path: str
    :param overview_outdir:
    :type overview_outdir:
    :param master_blastx_path: Path to the blastx output file
    :type master_blastx_path: str
    :return: annotated_clust_file: Path to the annotated mmseqs output file
    """

    mmseqs_outdir = os.path.join(overview_outdir, "mmseqs")
    if not os.path.exists(mmseqs_outdir):
        os.mkdir(mmseqs_outdir)

    mmseqs_outpath = run_mmseqs_cluster(master_contig_path, mmseqs_outdir)
    mmseqs_outdir = os.path.split(mmseqs_outpath)[0]

    for file in os.listdir(mmseqs_outdir):
        if file.endswith('cluster.tsv'):
            mmseqs_summary = file
            logger.log(f"mmseqs cluster file found: {file}")
        else:
            logger.log(f"Not the cluster file I am looking for: {file}")

    mmseqs_summary_path = os.path.join(mmseqs_outdir, mmseqs_summary)

    annotate_dir = os.path.join(overview_outdir, "annotate")
    if not os.path.exists(annotate_dir):
        os.mkdir(annotate_dir)

    if not os.path.exists(os.path.join(annotate_dir, f"master_contigs_mmseqs_annotated.tsv")):
        annotated_clust_file = blastx_annotate_mmseqs(mmseqs_summary_path, master_blastx_path, annotate_dir)
    else:
        annotated_clust_file = os.path.join(annotate_dir, f"master_contigs_mmseqs_annotated.tsv")

    filtered_annotations(annotated_clust_file, annotate_dir)

    return annotated_clust_file


def count_char_in_file(filename, char):
    """
    Counts the number of occurrences of a character in a file.

    :param filename: Path to the file
    :type filename: str
    :param char: Character to be counted
    :type char: str
    :return: count: Number of occurrences of the character
    :rtype: int
    """
    with open(filename, 'r') as file:
        content = file.read()
        count = 0
        for c in content:
            if c == char:
                count += 1
        return count


def parse_master_contigs(master_contig_path):
    """
    Parses the master contig file and returns a list of contig names.

    :param master_contig_path: Path to the master contig file
    :type master_contig_path: str
    :return: contig_list: List of contig names
    :rtype: list
    """

    contig_list = []
    with open(master_contig_path, 'r') as handle:
        for line in handle:
            if line.startswith('>'):
                contig = line.strip().split('>')[1]
                contig_list.append(contig)
    return contig_list


def get_meta_master_contigs(contigs_list):
    """
    Parses the master contig list and returns a dictionary containing the contig name, length and kmer coverage.

    :param contigs_list: List of contig names
    :type contigs_list: list
    :return: meta_master_dict: Dictionary containing the contig name, length and kmer coverage
    """

    meta_master_dict = {}

    for contig in contigs_list:
        tmp_line = contig.split('_')
        length = tmp_line[4]
        kmer_cov = tmp_line[6]

        if contig not in meta_master_dict:
            meta_master_dict[contig] = {'length': length, 'kmer_cov': kmer_cov}
        else:
            print("Duplicate contig found")

    return meta_master_dict


def write_metadata_file(master_contig_path, master_contig_path_trans, palm_all_file, mmseqs_outdir, blastx_set, hmm_set,
                        only_hmm_set):
    """
    Writes a metadata file containing the analysis results.

    :param master_contig_path: Path to the master contig file
    :type master_contig_path: str
    :param master_contig_path_trans:
    :type: str
    :param palm_all_file:
    :type: str
    :param mmseqs_outdir:
    :type: str
    :param blastx_set:
    :type: set
    :param hmm_set:
    :type: set
    :param only_hmm_set:
    :type: set
    :return: metadata_fn: Path to the metadata file
    :rtype: str
    """

    viral_contig_count = count_char_in_file(master_contig_path, ">")
    viral_contig_trans_count = count_char_in_file(master_contig_path_trans, ">")

    blastx_count = len(blastx_set)
    hmm_count = len(hmm_set)
    hmm_only = hmm_set.difference(blastx_set)
    blastx_only = blastx_set.difference(hmm_set)
    both = blastx_set.intersection(hmm_set)
    for file in os.listdir(mmseqs_outdir):
        if file.endswith('rep_seq.fasta'):
            mmseqs_rep_path = os.path.join(mmseqs_outdir, file)
            logger.log(f"mmseqs rep file found: {file}")

    mmseqs_rep_count = count_char_in_file(mmseqs_rep_path, ">")

    master_list = parse_master_contigs(master_contig_path)
    meta_master_dict = get_meta_master_contigs(master_list)

    length_500 = 0
    length_1000 = 0
    length_2000 = 0
    kmer_cov_5 = 0
    kmer_cov_10 = 0
    for contig, meta_dict in meta_master_dict.items():
        # TODO: Fix this
        if meta_dict['length'] != 'length':
            if float(meta_dict['length']) >= 500:
                length_500 += 1
            if float(meta_dict['length']) >= 1000:
                length_1000 += 1
            if float(meta_dict['length']) >= 2000:
                length_2000 += 1
            if float(meta_dict['kmer_cov']) >= 5:
                kmer_cov_5 += 1
            if float(meta_dict['kmer_cov']) >= 10:
                kmer_cov_10 += 1

    metadata_fn = f"{os.path.split(master_contig_path)[0]}/metadata.tsv"

    with open(metadata_fn, 'w') as out_handle:
        out_handle.write(f"Putative viral contigs (nucleotide): \t{viral_contig_count}\n")
        out_handle.write(f"Putative viral contigs (translated nucleotides): \t{viral_contig_trans_count}\n")
        out_handle.write(f"Contigs detected as viral only with HMMs: \t{len(only_hmm_set)}\n")
        out_handle.write(f"PalmScan (blastx only):\t{blastx_count}\n")
        out_handle.write(f"PalmScan (hmm only): \t{hmm_count}\n")
        out_handle.write(f"PalmScan sequences found in both: \t{len(both)}\n")
        out_handle.write(f"PalmScan sequences found in blastx only: \t{len(blastx_only)}\n")
        out_handle.write(f"PalmScan sequences found in hmm only: \t{len(hmm_only)}\n")
        out_handle.write(f"MMseqs2 clusters for HMM sequences, for .9 min-seq-id & c 0.8: \t{mmseqs_rep_count}\n")
        out_handle.write(f"Contigs >= 500 bp: \t{length_500}\n")
        out_handle.write(f"Contigs >= 1000 bp: \t{length_1000}\n")
        out_handle.write(f"Contigs >= 2000 bp: \t{length_2000}\n")
        out_handle.write(f"Contigs with kmer coverage >= 5: \t{kmer_cov_5}\n")
        out_handle.write(f"Contigs with kmer coverage >= 10: \t{kmer_cov_10}\n")

    return metadata_fn


def vir_disc_wrap(in_path):
    """
    Wrapper function for the viral discovery pipeline.

    :param in_path: Path to the input directory
    :type in_path: str
    :param blastx_dir_path: Path to the blastx output directory
    :type in_path: str
    :return: None
    """
    analysis_dict = {}


    samples = extract_dir_names(in_path)
    logger.log(f"Running viral discovery pipeline for samples: {' '.join(samples)}")

    out_dirs = create_outdir(in_path, samples)
    logger.log(f"Output directories created: {' '.join(out_dirs)}")

    general_outdir = os.path.split(out_dirs[0])[0]

    for sample, out_dir in zip(samples, out_dirs):

        in_dir = os.path.join(in_path, sample)
        contig_path = os.path.join(in_dir, 'contigs.fasta')

        # seqkit
        seqkit_outdir = os.path.join(out_dir, 'seqkit')
        if not os.path.exists(seqkit_outdir):
            os.mkdir(seqkit_outdir)
        logger.log(f"Running seqkit for sample {sample}")
        filter_contigs_path = os.path.join(os.path.join(seqkit_outdir, f"{sample}_seqkit.fasta"))
        if not os.path.exists(filter_contigs_path):
            run_seqkit_filt(filter_contigs_path, contig_path, sample)

        # transeq
        # Create transeq output directory
        logger.log(f"Creating transeq outdir {sample}")
        transeq_outdir = os.path.join(out_dir, 'transeq')
        if not os.path.exists(transeq_outdir):
            os.mkdir(transeq_outdir)
        logger.log(f"Running transeq for sample {sample}")
        transeq_out_path = os.path.join(transeq_outdir, f"{sample}_transeq.fasta")
        if not os.path.exists(transeq_out_path):
            run_transeq(transeq_out_path, filter_contigs_path, sample)

        # hmmsearch
        # Create hmm output directory
        hmm_outdir = os.path.join(out_dir, 'hmm')
        if not os.path.exists(hmm_outdir):
            os.mkdir(hmm_outdir)
        logger.log(f"Running hmmsearch for sample {sample} | rvmt & rdrpscan| 6-frame translation & gencode 1")
        analysis_dict = wrap_hmmsearch(transeq_out_path, hmm_rvmt, hmm_rdrp, hmm_outdir, sample, analysis_dict)


    # Summary
    overview_outdir = os.path.join(general_outdir, "overview_results")
    if not os.path.exists(overview_outdir):
        os.mkdir(overview_outdir)

    overview_disc_dir = os.path.join(overview_outdir, "viral_discovery")
    if not os.path.exists(overview_disc_dir):
        os.mkdir(overview_disc_dir)

    logger.log(f"Writing summary file...")
    sum_fn = write_summary_file(analysis_dict, overview_disc_dir)
    logger.log(f"Summary file written to {sum_fn}")
    # Filtered summary file
    logger.log(f"Writing filtered summary file...")
    filtered_sum_fn = write_filt_summary_file(sum_fn)
    logger.log(f"Filtered summary file written to {filtered_sum_fn}")

    # Extract contigs
    logger.log(f"Extracting contigs from filtered summary file...")
    master_contig_path = wrap_extract_contigs(in_path, overview_disc_dir, filtered_sum_fn)
    master_contig_path_trans = wrap_extract_trans(in_path, overview_disc_dir, general_outdir, filtered_sum_fn)
    logger.log(f"Contigs extracted to {master_contig_path}")

    # Highest scoring hmm hits
    logger.log(f"Writing the highest scoring hmm hits to file...")
    hmm_outfile = highest_scoring_hmms(filtered_sum_fn, overview_disc_dir)

    return overview_outdir,master_contig_path, hmm_outfile


