""" A radar chart (a.k.a. a spider or star chart) """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection


def radar_factory(num_vars, frame='circle'):
    """Create a radar chart with `num_vars` axes. This function creates a
    RadarAxes projection and registers it. """
    
    # calculate evenly-spaced axis angles
    theta = 2*np.pi * np.linspace(0, 1-1./num_vars, num_vars)
    # rotate theta such that the first axis is at the top
    theta += np.pi/2

    def draw_poly_patch(self):
        verts = unit_poly_verts(theta)
        return plt.Polygon(verts, closed=True, edgecolor='k')

    class RadarAxes(PolarAxes):

        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1
        # define draw_frame method
        draw_patch = draw_poly_patch

        def fill(self, *args, **kwargs):
            """Override fill so that line is closed by default"""
            closed = kwargs.pop('closed', True)
            return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super(RadarAxes, self).plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(theta * 180/np.pi, labels)

        def _gen_axes_patch(self):
            return self.draw_patch()

        def _gen_axes_spines(self):
            # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
            spine_type = 'circle'
            verts = unit_poly_verts(theta)
            # close off polygon by repeating first vertex
            verts.append(verts[0])
            path = Path(verts)

            spine = Spine(self, spine_type, path, visible = False)
            spine.set_transform(self.transAxes)
            return {'polar': spine}

    register_projection(RadarAxes)
    return theta


def unit_poly_verts(theta):
    """ Return vertices of polygon for subplot axes. This polygon is
    circumscribed by a unit circle centered at (0.5, 0.5) """
    
    x0, y0, r = [0.5] * 3
    verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    return verts



def plot(data, fig, colors, dormancy, dims, fs):
        
        """
        data  :  data for plot
        fig  :  matplotlib figure object
        dormancy   :  yes or no
        dims  :  
        fs  : fontsize """
        
        #if __name__ == '__main__':
        r, c, j = dims
        ax = fig.add_subplot(r, c, j)
        spoke_labels = data.pop('column names')
        N = len(spoke_labels)
        theta = radar_factory(N, frame='polygon')
        
        # Plot the data
        for keyname in data.keys():
            ax = fig.add_subplot(r, c, j, projection='radar')
            plt.rgrids([1])
            for d, color in zip(data[keyname], colors):
                ax.plot(theta, d, color=color, alpha=0.6)
                ax.fill(theta, d, facecolor=color, alpha=0.2)
            ax.set_varlabels(spoke_labels)
        
        if dormancy == False:
            plt.figtext(0.5, 0.965, '8-Factor Hydrobide model w/out dormancy',
                    ha='center', color='black', weight='bold', size='large')
        else:
            plt.figtext(0.5, 0.965, '8-Factor Hydrobide model with dormancy',
                    ha='center', color='black', weight='bold', size='large')
        
        return fig