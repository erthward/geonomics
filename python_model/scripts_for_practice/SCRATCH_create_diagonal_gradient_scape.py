from scipy import interpolate



def create_gradient(dim):
    dim_plus = dim + 1
    dub_dim_plus = dim_plus*2
    x = list(range(dim_plus)) + [dim]*dim_plus + list(reversed(list(range(dim_plus)))) + [0] *dim_plus
    y = [dim]*dim_plus + list(reversed(list(range(dim_plus)))) + [0] * dim_plus + list(range(dim_plus))
    vals = list(range(dim_plus)) + list(range(dim_plus, dub_dim_plus)) + list(reversed(list(range(dim_plus,dub_dim_plus)))) + list(reversed(list(range(dim_plus))))
    vals = [i/max(vals) for i in vals]
    pts = np.vstack((np.array(x), np.array(y))).T
    gx, gy = np.mgrid[1:dim:complex("%ij"%dim), 1:dim:complex("%ij"%dim)]
    out = interpolate.griddata(pts, vals, (gx, gy), method = 'cubic')
    return(out)
