#!/usr/bin/env python
# Should be python3 in the near future
from __future__ import division
from pylab import *
fig, ax = subplots(figsize=(6, 1))
fig.subplots_adjust(bottom=0.5)


colors = [array([255,  140,  255]) / 255.,  array([221,  111,  242]) / 255., array([187,  82,  229]) / 255.,  array([153,  53,  216]) / 255.,  array([119,  24,  203]) / 255.,  array([0,  0,  187]) / 255.,  array([0,  44,  204]) / 255.,  array([0,  88,  221]) / 255.,  array([0,  132,  238]) / 255.,  array([0,  175,  255]) / 255.,  array([0,  235,  255]) / 255.,  array([39,  255,  215]) / 255.,  array([99,  255,  155]) / 255.,  array([163,  255,  91]) / 255.,  array([211,  255,  43]) / 255.,  array([255,  255,  0]) / 255.,  array([255,  207,  0]) / 255.,  array([255,  159,  0]) / 255.,  array([255,  111,  0]) / 255.,  array([255,  63,  0]) / 255.,  array([255,  0,  0]) / 255.,  array([216,  0,  15]) / 255.,  array([178,  0,  31]) / 255.,  array([140,  0,  47]) / 255.,  array([102,  0,  63]) / 255.,  array([52,  52,  52]) / 255.,  array([96,  96,  96]) / 255.,  array([140,  140,  140]) / 255.,  array([184,  184,  184]) / 255.,  array([228,  228,  228]) / 255., [1.,1.,1.] ]


cmap = mpl.colors.ListedColormap(colors)
		
#		['red', 'green', 'blue', 'cyan'])
#cmap.set_under(array([255,255,255]) / 255.)
cmap.set_under([1,1,1])
cmap.set_over([0,0,0])

bounds = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 125, 150, 200, 300,600]

print len(bounds), len(colors)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb2 = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, boundaries=[0] + bounds + [600], extend='both', ticks=bounds, spacing='proportional', orientation='horizontal')
cb2.set_label("O3 [PPBV]") #( Discrete intervals, some other units')
show()
