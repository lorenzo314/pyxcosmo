import matplotlib.pyplot as plt

def plotImageDiff(map_python, map_idl, cmap, label, save=False, scale=1.0, outdir='figures/'):
   diff = map_idl - map_python
   diff.shape, diff.min(), diff.max()
   diff *= scale

   fig=plt.figure()

   a=fig.add_subplot(1,3,1)
   imgplot = plt.imshow(map_python, origin='lower', cmap=cmap)
   a.set_title(label + ' python')
   plt.colorbar(orientation ='horizontal', ticks=[map_python.min(), map_python.max()])

   a=fig.add_subplot(1,3,2)
   imgplot = plt.imshow(map_idl, origin='lower', cmap=cmap)
   a.set_title(label + ' IDL')
   plt.colorbar(orientation ='horizontal', ticks=[map_idl.min(), map_idl.max()])

   a=fig.add_subplot(1,3,3)
   imgplot = plt.imshow(diff, origin='lower', cmap=cmap)
   plt.colorbar(orientation ='horizontal', ticks=[diff.min(), diff.max()])
   a.set_title("%s %.1e" % (label+'IDL-python', scale))

   if save: plt.savefig(outdir + 'map_' + label + '_idl_py_color32.png')

   return
