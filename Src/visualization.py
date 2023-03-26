import plotly.graph_objects as go

def read_material(filename):
    file = open(filename, 'r')
    
    prds = [float(prd_str) for prd_str in file.readline().split()]
    dims = [int(dim_str) for dim_str in file.readline().split()]
    
    incr = [prds[idx] / (dims[idx] - 1) for idx in range(3)]
    
    x = []
    y = []
    z = []
    values = []
    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                x.append(incr[0]*i)
                y.append(incr[1]*j)
                z.append(incr[2]*k)
                values.append(float(file.readline()))
            
    file.close()
    
    return x, y, z, values


def visualize(x, y, z, values, 
              x_show=True, y_show=True, z_show=True, 
              surface_count=2, 
              x_slices=[], y_slices=[], z_slices=[], 
              opacity=1,
              renderer="browser"):
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

    fig.show(renderer=renderer)















