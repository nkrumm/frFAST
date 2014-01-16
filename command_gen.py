import argparse
import pandas
import sys
import os

TEMPLATE_HEADER="""
#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load modules modules-init modules-gs modules-eichler
module load python/2.7.2
module load zlib/latest
module load hdf5/1.8.4
module load numpy/1.6.1
module load zeromq/2.1.11
"""

# TEMPLATE_ROW="""
# python {SCRIPT} \
# --sampleID={SAMPLE_ID} \
# --source={SOURCE_BAM} \
# --output={OUTPUT} \
# --log_dir={LOG_DIR} \
# --index_dir={INDEX_DIR} \
# --index={INDEX_NAME} \
# --translate_table={TRANSLATE_TABLE} \
# --port {PORT} --timeout 1800 \
# {OPTIONS} \
# --error_log={FR_FAST_LOG} 2>> {FR_FAST_LOG}
# """

TEMPLATE_ROW="""
if [ ! -e {OUTPUT} ]; then\n \
python {SCRIPT} \
--sampleID={SAMPLE_ID} \
--source={SOURCE_BAM} \
--output={OUTPUT} \
--log_dir={LOG_DIR} \
--index_dir={INDEX_DIR} \
--index={INDEX_NAME} \
--translate_table={TRANSLATE_TABLE} \
--port {PORT} --timeout 1800 \
{OPTIONS} \
--error_log={FR_FAST_LOG}
fi
"""

FRFAST_SCRIPT="/net/eichler/vol8/home/nkrumm/EXOMES/frFAST/controller2.1.py"

if __name__ == '__main__':
    parser = argparse.ArgumentParser("""This script can generate batch files for processing large numbers of exomes with frFAST.
        The input is a tab-delimited text file with the following columns:

        sampleID    bam_path
        sample_1    /path/to/bam/file/sample_1.bam
        sample_B    /path/to/bam/file/sample_B.bam
        ...

        Example usage for creating a batch file suitable for running a headless (--disable-gui) frFAST on a single node (--single-host):

        python command_gen.py \\
            sample_table.txt \\
            /path/to/output_dir/ \\
            /path/to/log_dir/ \\
            /net/grc/shared/scratch/nkrumm/INDEX/default_exome/default_exome.fa \\
            /net/grc/shared/scratch/nkrumm/translate_tables/default_exome.translate.txt \\
            --dont-rsync-index --disable-gui --single-host --disable-port-scan 

        Please note that you may have to adjust the "module load" and SGE directives in the TEMPLATE_HEADER variable.
        The defaults are:

        #!/bin/bash
        #$ -S /bin/bash
        #$ -cwd

        module load modules modules-init modules-gs modules-eichler
        module load python/2.7.2
        module load zlib/latest
        module load hdf5/1.8.4
        module load numpy/1.6.1
        module load zeromq/2.1.11

    """)
    parser.add_argument("input_table", help="Tab-delimited table with header field containing columns [sampleID, bam_path]")
    parser.add_argument("output_directory", help="Path for output frFAST files")
    parser.add_argument("log_directory", help="Path for frFAST log files")
    parser.add_argument("index_fasta", help="Path in to mrsFAST compatible index.fa file for frFAST to use")
    parser.add_argument("translate_table", help="Path in which frFAST will look for a translate_table file")
    parser.add_argument("--port", default=8000, help="Start of port range for frFAST to use")
    parser.add_argument("--disable-port-scan", action="store_true", default=False, help="See frFAST documentation")
    parser.add_argument("--single-host",  action="store_true", default=False, help="See frFAST documentation")
    parser.add_argument("--dont-rsync-index",  action="store_true", default=False, help="See frFAST documentation")
    parser.add_argument("--disable-gui",  action="store_true", default=False, help="See frFAST documentation")
    parser.add_argument("-o", "--output-file", help="Batch output file (default is stdout)", default="STDOUT")
    args = parser.parse_args()
    samples = pandas.read_csv(args.input_table, sep="\t")
    if ("sampleID" not in samples) or ("bam_path" not in samples):
        raise("Error, input_table must contain 'sampleID' and 'bam_path' fields and must be tab-delimited!")

    if args.output_file != "STDOUT":
        out_f = open(args.output_file, 'w')
    else:
        out_f = sys.stdout

    out_f.write(TEMPLATE_HEADER)
    

    for ix, row in samples.iterrows():
        if "output_path" not in row:
            output_path = os.path.join(args.output_directory, "%s.h5" % row["sampleID"])
        else:
            output_path = row["output_path"]

        log_dir = os.path.join(args.log_directory, row["sampleID"])
        
        options = filter(lambda opt: vars(args)[opt] == True, ["disable_port_scan", "single_host", "dont_rsync_index", "disable_gui"])
        opt_string = ""
        for opt in options:
            opt_string += "--%s " % opt.replace("_", "-")

        opt_string = opt_string.lstrip(" ")

        out_f.write(TEMPLATE_ROW.format(
            SCRIPT=FRFAST_SCRIPT,
            SOURCE_BAM=row["bam_path"],
            OUTPUT=output_path,
            LOG_DIR=log_dir,
            SAMPLE_ID=row["sampleID"],
            INDEX_DIR=os.path.dirname(args.index_fasta),
            INDEX_NAME=os.path.basename(args.index_fasta),
            TRANSLATE_TABLE=args.translate_table,
            PORT=args.port,
            OPTIONS=opt_string,
            FR_FAST_LOG=os.path.join(log_dir, "frfast.errors.log")
            ))

    out_f.close()