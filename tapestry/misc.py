from ._version import __version__

import os, sys, argparse, itertools, errno
import logging as log
from functools import partial, lru_cache

from plumbum import local, CommandNotFound

failed = []
tools = {'paftools':'paftools.js',
         'minimap2':'minimap2', 
         'mosdepth':'mosdepth', 
         'samtools':'samtools', 
         'zgrep':'zgrep', 
         'pigz':'pigz',
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


class PAF:
    def __init__(self, pafline):
        f = pafline.rstrip().split('\t')
        self.query_name      = f[0]
        self.query_length    = int(f[1])
        self.query_start     = int(f[2])
        self.query_end       = int(f[3])
        self.strand          = f[4]
        self.subject_name     = f[5]
        self.subject_length   = int(f[6])
        self.subject_start    = int(f[7])
        self.subject_end      = int(f[8])
        self.matches         = int(f[9])
        self.block_length    = int(f[10])
        self.mapping_quality = int(f[11])


def get_args(arglist=[], description="", scriptargs=[]):

    parser = argparse.ArgumentParser(description=description)
    for scriptarg in scriptargs:
        exec(f'parser.add_argument({scriptarg})')

    parser.add_argument('-c', '--cores', help="number of parallel cores to use", type=int, default=1)
    parser.add_argument('-v', '--verbose', help="report on progress", action="count", default=0)
    parser.add_argument('-V', '--version', help="report version number and exit", action="store_true")

    # If no arguments, print usage message
    if not arglist:
        parser.print_help()
        sys.exit()

    args = parser.parse_args(arglist)

    if args.version:
        versions()
        sys.exit()

    log.basicConfig(format="%(asctime)s %(levelname)s\t%(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    set_verbosity(args.verbose)

    if args.cores < 1:
        log.error("Must specify at least one core")
        sys.exit()

    cores_available = len(os.sched_getaffinity(0)) # https://docs.python.org/3/library/os.html#os.cpu_count
    if args.cores > cores_available:
        log.error(f"{args.cores} cores requested but only {cores_available} available, please reduce cores")
        sys.exit()

    return args


def get_weave_args(arglist=[]):
    args = get_args(sys.argv[1:], 
           "weave: assess quality of one genome assembly",
           ["'-a', '--assembly', help='filename of assembly in FASTA format', type=str",
            "'-r', '--reads', help='filename of reads in FASTQ format (can be gzipped)', type=str",
            "'-t', '--telomere', help='telomere sequence to search for', type=str, action='append', nargs='+'",
            "'-o', '--output', help='directory to write output, default weave_output', type=str, default='weave_output'"])

    if not args.assembly:
        log.error("Assembly file in FASTA format is required (-a, --assembly)")
        sys.exit()

    return args


def get_stitch_args(arglist=[]):
    args = get_args(sys.argv[1:], # Pass in arguments so this function can be unit tested
           "stitch: compare genome assemblies",
           ["'-a', '--assemblies', help='tapestry folders', type=str, action='append', nargs='+'",
            "'-o', '--output', help='directory to write output, default stitch_output', type=str, default='stitch_output'"]) 

    if len(args.assemblies[0]) < 2:
        log.error("Please specify at least two Tapestry folders to process (-a, --assemblies)")
        sys.exit()

    return args


def setup_output(outdir):
    try:
        os.mkdir(outdir)
        log.info(f"Created output directory {outdir}")
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            log.warning(f"Output directory {outdir} found, will use existing analysis files if present, but overwrite reports")
        else:
            raise



def versions(verbosity=2):

    log.getLogger().setLevel(log.INFO) # Suppress plumbum DEBUG messages
    print(f"Tapestry version {__version__}")

    debug = "Dependencies\n"

    if 'minimap2' in globals():
        debug += f"minimap2\t{minimap2('--version').rstrip()}\t{minimap2}\n"
    else:
        debug += f"minimap2\tMISSING\n"

    if 'paftools' in globals():
        debug += f"paftools\t{paftools('version').rstrip()}\t{paftools}\n"
    else:
        debug += f"paftools\tMISSING\n"

    samtools_version = samtools['--version'] | head['-n 1'] | cut['-d ', '-f2']
    debug += f"samtools\t{samtools_version().rstrip()}\t{samtools}\n"

    mosdepth_version = mosdepth['-h'] | head['-n 1'] | cut['-d ', '-f2']
    debug += f"mosdepth\t{mosdepth_version().rstrip()}\t{mosdepth}\n"

    set_verbosity(verbosity) # Reset logger now plumbum commands are done

    print(debug.expandtabs(15))


def set_verbosity(verbosity):
    if verbosity == 1:
        log.getLogger().setLevel(log.INFO)
    elif verbosity > 1:
        log.getLogger().setLevel(log.DEBUG)
    else:
        log.getLogger().setLevel(log.WARN)


def flatten(l):
    return list(itertools.chain.from_iterable(l))


def grep(pattern, filename):
    return zgrep(pattern, filename, retcode=(0,1)).split('\n')[:-1]


def cached_property(function):
    return property(lru_cache()(function))
