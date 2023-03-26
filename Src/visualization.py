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


def visualize(x, y, z, values):
    return




