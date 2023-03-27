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


def main(filename_in, filename_out):
    x, y, z, values = read_material(filename_in)
    
    # specify visualization parameters
    visualize(x, y, z, values, filename=filename_out)





main(sys.argv[1], sys.argv[2])


