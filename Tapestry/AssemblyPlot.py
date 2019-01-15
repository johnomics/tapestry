import logging as log

from plotly.offline import plot
from plotly.graph_objs import Bar, Layout, Scatter

from .misc import flatten

class AssemblyPlot:

    def _assemblyplot(self, name, plotobj):
        plotfile = f"{self.outdir}/{name}.html"
        log.debug(f"Making plot {plotfile}")
        try:
            plot(plotobj, auto_open=False, filename=plotfile)
        except:
            log.error(f"Can't make plot {plotfile}")
            sys.exit()

    def lengthplot(self):
        self._assemblyplot("lengths", {
                "data": [Bar(x=[c.name for c in self.contigs], y=[len(c) for c in self.contigs])],
                "layout": Layout(title="Contig Lengths")
        })

    def depthplot_full(self):
        self._assemblyplot("read_depths_full", {
                "data": [Scatter(x=[d.end for d in c.depths('reads')], y=[d.depth for d in c.depths('reads')], name=c.name, mode='lines') for c in self.contigs],
                "layout": Layout(title="Contig Read Depths")
        })

    def depthplot_summary(self):
        xnames     = flatten([[c.name]    * len(c.depths('reads')) for c in self.contigs])
        ydepths    = flatten([[d.depth for d in c.depths('reads')] for c in self.contigs])

        self._assemblyplot("read_depths_summary", {
                "data": [Scatter(x=xnames, y=ydepths, mode="markers", marker=dict(opacity=0.1))],
                "layout": Layout(title="Contig Read Depth Summary", shapes=[
                    {'type':'line', 'x0':xnames[0], 'x1':xnames[-1], 'y0':self.median_depth, 'y1':self.median_depth},
                    {'type':'line', 'x0':xnames[0], 'x1':xnames[-1], 'y0':self.median_depth*2, 'y1':self.median_depth*2},
                    {'type':'line', 'x0':xnames[0], 'x1':xnames[-1], 'y0':self.median_depth/2, 'y1':self.median_depth/2}
                ]) 
        })
