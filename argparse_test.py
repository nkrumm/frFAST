import argparse
parser = argparse.ArgumentParser(description='This is frFAST (frankenFAST), a lightweight pipeline designed to calculate read-depth across an exome or genome using the mrsFAST mapper. \nNik Krumm, 2011')
parser.add_argument('--source','-s', metavar='/path/to/source_file.bam', type=str, nargs='?',help="Source BAM to be processed")
parser.add_argument('--output','-o', metavar='/path/to/output_file.h5', type=str, nargs='?',help="Output location of HDF5 file")
parser.add_argument('--port','-p', metavar='5555', type=int, nargs='?', default = 5555,\
	help="TCP port offset to use. Will use FOUR CONSECUTIVE ports starting with this value. Default port = 5555")
parser.add_argument('--disable_port_scan',action='store_true',\
	help="Disable initial port scan. Use with caution and only with ports you know are open!")

args = parser.parse_args()
print args.port
print args.disable_port_scan
