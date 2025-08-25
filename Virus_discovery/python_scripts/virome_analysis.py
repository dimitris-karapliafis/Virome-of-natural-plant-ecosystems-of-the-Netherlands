from vir_disc_module import *
from annotation_module import *
from filtering_module import *
from logger import logger


##########        DATABASES     ##########

hmm_rdrp = "PHD/DBs/RDRP-scan/RdRp_HMM_profile_CLUSTALO.db"
hmm_neordrp = "PHD/DBs/neo-rdrp/NeoRdRp-HMM.v1.1.hmm"
hmm_firth = "/DBs/firth/conc_prof.hmm"
hmm_rvmt = "PHD/DBs/rvmt/RVMT.hmm"
diam_rdrpscan = "PHD/DBs/RDRP-scan/RdRp-scan_0.90.dmnd"
diam_neordrp = "PHD/DBs/neo-rdrp/NeoRdRp-seq.v1.1.dmnd"
genomad_db = "PHD/DBs/genomad_db"
dblast_nr = 'PHD/DBs/blast_nr/diamond_nr.dmnd'


##########        ARGUMENTS     ##########
parser = argparse.ArgumentParser(description='VirDisc pipeline', usage= 'vir_disc_pipeline.py -i <input_dir>')
parser.add_argument('-i', '--input', help='Path to the input directory', required=True)
args = parser.parse_args()


if __name__ == "__main__":
    # Create logger
    logger = Logger('vir_disc_pipeline.log')
    logger.log("Starting the script...")
    logger.start_timer()

    # Parse arguments
    args = parser.parse_args()

    # Run the virus discovery module
    overview_outdir, master_contigs_path,high_score_hmm_fn = vir_disc_wrap(args.input)
    logger.log("Virus discovery module finished.")
    logger.log("Starting the annotation module...")
    mmseqs_outpath_clust,blast_outpath_plantvir, palmscan_outpath,blast_outpath_palmdb, blast_outpath_nr, mmseqs_outpath_embryo, genomad_outpath, mmseqs_nt_fn= wrap_annotation(master_contigs_path, overview_outdir)
    annotated_fn = wrap_annotation_file(overview_outdir,mmseqs_outpath_clust, blast_outpath_plantvir, palmscan_outpath, blast_outpath_palmdb,
                    blast_outpath_nr, mmseqs_outpath_embryo, genomad_outpath, high_score_hmm_fn,mmseqs_nt_fn)

    filtered_fn = wrap_filtering(overview_outdir, annotated_fn, mmseqs_nt_fn,master_contigs_path)




    # Stop the logger
    logger.stop_timer()
    logger.log("Script finished.")
