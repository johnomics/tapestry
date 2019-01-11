import logging as log

from ._version import __version__

from .misc import set_verbosity

from plumbum import local
from plumbum.cmd import minimap2, mosdepth, samtools
from plumbum.cmd import head, cut, zgrep, pigz
paftools = local["paftools.js"]

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
