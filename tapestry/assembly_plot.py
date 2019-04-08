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
                "data": [Bar(x=[self.contigs[c].name for c in self.contiglist], 
                             y=[len(self.contigs[c]) for c in self.contiglist])],
                "layout": Layout(title="Contig Lengths")
        })

    def depthplot_full(self):
        self._assemblyplot("read_depths_full", {
                "data": [Scatter(x=self.contigs[c].read_depths['end'],
                                 y=self.contigs[c].read_depths['depth'], 
                                 name=self.contigs[c].name, mode='lines') for c in self.contiglist],
                "layout": Layout(title="Contig Read Depths")
        })

    def depthplot_summary(self):
        xnames     = flatten([[self.contigs[c].name]    * len(self.contigs[c].read_depths.index) for c in self.contiglist])
        ydepths    = flatten([self.contigs[c].read_depths['depth'] for c in self.contiglist])

        self._assemblyplot("read_depths_summary", {
                "data": [Scatter(x=xnames, y=ydepths, mode="markers", marker=dict(opacity=0.1))],
                "layout": Layout(title="Contig Read Depth Summary", shapes=[
                    {'type':'line', 'x0':xnames[0], 'x1':xnames[-1], 'y0':self.median_depth, 'y1':self.median_depth},
                    {'type':'line', 'x0':xnames[0], 'x1':xnames[-1], 'y0':self.median_depth*2, 'y1':self.median_depth*2},
                    {'type':'line', 'x0':xnames[0], 'x1':xnames[-1], 'y0':self.median_depth/2, 'y1':self.median_depth/2}
                ]) 
        })
