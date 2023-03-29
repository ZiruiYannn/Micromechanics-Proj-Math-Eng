import plotly.graph_objects as go
import sys

# read material information from file
def read_material(filename):
    file = open(filename, 'r')
    
    # read in periods which corresponds to the physical length of each dimension
    prds = [float(prd_str) for prd_str in file.readline().split()]
    # read the number of voxels in each dimension 
    dims = [int(dim_str) for dim_str in file.readline().split()]
    
    incr = [prds[idx] / (dims[idx] - 1) for idx in range(3)]
    
    x = []
    y = []
    z = []
    values = []
    # read in value for every voxel and process data for visualization 
    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                x.append(incr[0]*i)
                y.append(incr[1]*j)
                z.append(incr[2]*k)
                values.append(float(file.readline()))
            
    file.close()
    
    return x, y, z, values


# visualize material for given discretization and values
# arguments:
#           x, y, z: lists containing the coordinates of each value in their respective dimension
#           values: list of values at each voxel corresponding to coordinates
#           x_show, y_show, z_show: booleans specifying whether to visualize outside surface(s) for respective dimension
#           surface_count: number of equidistant material slices to visualize in every dimension
#           x_slices, y_slices, z_slices: array indicating the points at which to visualize a material slice for respective dimension 
#                                         empty array results in no specific material slice to be visualized for respective dimension 
#           opacity: opacity level of visualization with 1 being non-opaque and 0 being completely seetrough
#           filename: file location and name + extension to save the visulization on disk
#                     None results in no image locally stored
#           renderer: plotly renderer to use when visualizing the material (browser, png, vsg, ...)
#                     None results in no visualization being displayed
def visualize(x, y, z, values, 
              x_show=True, y_show=True, z_show=True, 
              surface_count=2, 
              x_slices=[], y_slices=[], z_slices=[], 
              opacity=1,
              filename=None,
              renderer=None):
    if len(x_slices) == 0: 
        x_slices_show = False
    if len(y_slices) == 0: 
        y_slices_show = False
    if len(z_slices) == 0: 
        z_slices_show = False
    
    fig= go.Figure(data=go.Isosurface(
            x = x,
            y = y,
            z = z,
            value = values,
            caps = dict(x_show = x_show, y_show = y_show, z_show = z_show),
            surface_count = surface_count,
            slices_x = dict(show=x_slices_show, locations=x_slices),
            slices_y = dict(show=y_slices_show, locations=y_slices),
            slices_z = dict(show=z_slices_show, locations=z_slices),
            opacity = opacity,
        ))
    
    if filename != None:
        fig.write_image(filename)
    
    if renderer != None:
        fig.show(renderer=renderer)





# main
# input: input file, output file, [options]
# visualization options:
#                       -disp: displays interactive visualization in browser
#                              this option can be combined with the other options, however, if specified, it should always be the first 
#                       -surf: displays a numberof equidistant surfaces in a given dimensio
#                              parameters:
#                                          dimension: x, y, or z
#                                          surface_count: number of equidistant surfaces to visualize
#                                          opacity: the opacity of the visualization
#                                                   real number between 0 and 1, with 1 being non-opaque and 0 being completely seetrough 
#                       -slice: displays slices in a specified dimension
#                               parameters:
#                                          dimension: x, y or z
#                                          slices: a variable number of real numbers specifying the position of the slice in the specified dimension
#                                                   should be in the range of the read in material lengths
#                                          opacity: the opacity of the visualization
#                                                   real number between 0 and 1, with 1 being non-opaque and 0 being completely seetrough 
filename_in = sys.argv[1]
filename_out = sys.argv[2]
x, y, z, values = read_material(filename_in)

pad = 0
renderer = None
if len(sys.argv) > 3:
    if sys.argv[3] == "-disp":
        renderer = 'browser'
        
        if len(sys.argv)  > 4:
            pad = 1
        
    if sys.argv[pad+3] == "-surf":
        surf_count = sys.argv[pad+5]
        opct = sys.argv[pad+6]
        
        if sys.argv[pad+4] == "x":
            visualize(x, y, z, values, x_show=True, surface_count=surf_count, opacity=opct, filename=filename_out, renderer=renderer)
        
        if sys.argv[pad+4] == "y":
            visualize(x, y, z, values, y_show=True, surface_count=surf_count, opacity=opct, filename=filename_out, renderer=renderer)
        
        if sys.argv[pad+4] == "z":
            visualize(x, y, z, values, z_show=True, surface_count=surf_count, opacity=opct, filename=filename_out, renderer=renderer) 
            
    elif sys.argv[pad+3] == "-slice":
        opct = sys.argv[-1]
        
        if sys.argv[pad+4] == "x":
            slices = []
            for i in range(pad+5, len(sys.argv)-1):
                slices.append(sys.argv[pad+5+i])
            
            visualize(x, y, z, values, x_slices=slices, opacity=opct, filename=filename_out, renderer=renderer)
        
        if sys.argv[pad+4] == "y":
            slices = []
            for i in range(pad+5, len(sys.argv)-1):
                slices.append(sys.argv[pad+5+i])
            
            visualize(x, y, z, values, y_slices=slices, opacity=opct, filename=filename_out, renderer=renderer)
        
        if sys.argv[pad+4] == "z":
            slices = []
            for i in range(pad+5, len(sys.argv)-1):
                slices.append(sys.argv[pad+5+i])
            
            visualize(x, y, z, values, z_slices=slices, opacity=opct, filename=filename_out, renderer=renderer)
            
    else:
        print("unsupported command given")


















