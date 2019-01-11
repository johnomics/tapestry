import sys
import logging as log
from ._version import __version__
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
