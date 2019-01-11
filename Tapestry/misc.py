from ._version import __version__

import os, sys, argparse, itertools
import logging as log
from functools import partial

from plumbum import local, CommandNotFound

failed = []
tools = {'paftools':'paftools.js',
         'minimap2':'minimap2', 
         'mosdepth':'mosdepth', 
         'samtools':'samtools', 
         'zgrep':'zgrep', 
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



def get_args(arglist=[]):

    parser = argparse.ArgumentParser(description="Tapestry: assess genome assembly quality")

    parser.add_argument('-a', '--assembly', help="filename of assembly in FASTA format", type=str)
    parser.add_argument('-r', '--reads', help="filename of reads in FASTQ format (can be gzipped)", type=str)
    parser.add_argument('-t', '--telomere', help="telomere sequence to search for", type=str, action='append', nargs='+')
    parser.add_argument('-c', '--cores', help="number of parallel cores to use", type=int, default=1)
    parser.add_argument('-o', '--output', help="directory to write output, default tapestry_output", type=str, default="tapestry_output")
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

    log.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

    set_verbosity(args.verbose)

    if args.cores < 1:
        log.error("Must specify at least one core")
        sys.exit()

    cores_available = len(os.sched_getaffinity(0)) # https://docs.python.org/3/library/os.html#os.cpu_count
    if args.cores > cores_available:
        log.error(f"{args.cores} cores requested but only {cores_available} available, please reduce cores")
        sys.exit()

    if not args.assembly:
        log.error("Assembly file in FASTA format is required (-a, --assembly)")
        sys.exit()


    return args



def versions(verbosity=2):

    log.getLogger().setLevel(log.INFO) # Suppress plumbum DEBUG messages
    print(f"Tapestry version {__version__}")

    debug = "Dependencies\n"

    debug += f"minimap2\t{minimap2('--version').rstrip()}\t{minimap2}\n"

    debug += f"paftools\t{paftools('version').rstrip()}\t{paftools}\n"

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



#http://code.activestate.com/recipes/577452-a-memoize-decorator-for-instance-methods/
class memoize(object):
    """cache the return value of a method
    
    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments passed to a method decorated with memoize must
    be hashable.
    
    If a memoized method is invoked directly on its class the result will not
    be cached. Instead the method will be invoked like a static method:
    class Obj(object):
        @memoize
        def add_to(self, arg):
            return self + arg
    Obj.add_to(1) # not enough arguments
    Obj.add_to(1, 2) # returns 3, result is not cached
    """
    def __init__(self, func):
        self.func = func
    def __get__(self, obj, objtype=None):
        if obj is None:
            return self.func
        return partial(self, obj)
    def __call__(self, *args, **kw):
        obj = args[0]
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        key = (self.func, args[1:], frozenset(kw.items()))
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kw)
        return res
