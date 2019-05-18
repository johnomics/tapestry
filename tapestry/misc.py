from ._version import __version__

import os, sys, argparse, itertools, errno, io, pkg_resources
import logging as log
from functools import partial, lru_cache
from tqdm import tqdm

from plumbum import local, CommandNotFound

failed = []
tools = {'minimap2':'minimap2', 
         'samtools':'samtools', 
         'head':'head', 
         'cut':'cut'}


for tool in tools:
    try:
        exec(f"{tool} = local['{tools[tool]}']")
    except CommandNotFound:
        failed.append(tools[tool])


if failed:
    dep = 'dependency' if len(failed)==1 else 'dependencies'
    itthey = 'it is' if len(failed)==1 else 'they are'
    print(f"Tapestry can't find the following {dep}. Please check {itthey} installed and try again.")
    print('\n'.join(failed))
    sys.exit()


report_folder = pkg_resources.resource_filename(__name__, 'report')

tapestry_tqdm = partial(tqdm, unit=" contig", leave=False, miniters=1, dynamic_ncols=True)

def get_args(arglist=[], description="", scriptargs=[]):

    parser = argparse.ArgumentParser(description=description)
    for scriptarg in scriptargs:
        exec(f'parser.add_argument({scriptarg})')

    parser.add_argument('-c', '--cores', help="number of parallel cores to use (default 1)", type=int, default=1)
    parser.add_argument('-v', '--version', help="report version number and exit", action="store_true")

    # If no arguments, print usage message
    if not arglist:
        parser.print_help()
        sys.exit()

    args = parser.parse_args(arglist)

    if args.version:
        versions()
        sys.exit()

    log.basicConfig(format="%(asctime)s %(levelname)s\t%(message)s", datefmt="%Y-%m-%d %H:%M:%S", level=log.INFO)

    if args.cores < 1:
        log.error("Please specify at least one core")
        sys.exit()

    return args


def get_weave_args(arglist=[]):
    args = get_args(sys.argv[1:], 
           "weave: assess quality of one genome assembly",
           ["'-a', '--assembly', help='filename of assembly in FASTA format', type=str",
            "'-r', '--reads', help='filename of reads in FASTQ format (can be gzipped)', type=str",
            "'-d', '--depth', help='genome coverage to subsample from FASTQ file (default 50)', type=int, default=50",
            "'-l', '--length', help='minimum read length to retain when subsampling (default 10000)', type=int, default=10000",
            "'-t', '--telomere', help='telomere sequence to search for', type=str, action='append', nargs='+'",
            "'-w', '--windowsize', help='window size for ploidy calculations (default 10000)', type=int, default=10000",
            "'-o', '--output', help='directory to write output, default weave_output', type=str, default='weave_output'"])

    if not args.assembly:
        log.error("Assembly file in FASTA format is required (-a, --assembly)")
        sys.exit()

    return args


def weave_welcome(arglist=[]):
    versions()
    
    print("\nWelcome to Tapestry!\n")
    print(f"Assembly to validate\t{arglist.assembly}")
    print(f"Reads to sample from\t{arglist.reads}")
    print(f"Coverage to sample\t{arglist.depth}")
    print(f"Minimum read length\t{arglist.length}")
    print(f"Telomere sequence(s)\t{' '.join(arglist.telomere[0])}")
    print(f"Ploidy window size\t{arglist.windowsize}")
    print(f"Output directory\t{arglist.output}")
    print()


def versions():

    log.getLogger().setLevel(log.INFO) # Suppress plumbum DEBUG messages
    print(f"Tapestry version {__version__}")

    version_message = "Dependencies\n"

    if 'minimap2' in globals():
        version_message += f"minimap2\t{minimap2('--version').rstrip()}\t{minimap2}\n"
    else:
        version_message += f"minimap2\tMISSING\n"

    samtools_version = samtools['--version'] | head['-n 1'] | cut['-d ', '-f2']
    version_message += f"samtools\t{samtools_version().rstrip()}\t{samtools}\n"

    print(version_message.expandtabs(15))


def setup_output(outdir):
    try:
        os.mkdir(outdir)
        log.info(f"Created output directory {outdir}")
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            log.warning(f"Output directory {outdir} found, will use existing analysis files if present and up-to-date, but overwrite reports")
        else:
            raise


def file_exists(filename, deps=[]):
    return (os.path.exists(filename) and 
             all([os.stat(filename).st_mtime > os.stat(dep).st_mtime for dep in deps])
           )


def cached_property(function):
    return property(lru_cache()(function))


def include_file(filename):
    try:
        with io.open(os.path.join(report_folder, filename), "r", encoding='utf-8') as f:
            return f.read()
    except (OSError, IOError) as e:
        log.error(f"Could not include file '{filename}' in report: {e}")