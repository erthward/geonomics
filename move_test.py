import geonomics as gnx
mod = gnx.make_model('test_params_no_movement.py')
mod.walk(mode = 'burn', T = 1e6, verbose = True)
positions = {i:(v.x, v.y) for i, v in mod.comm[0].items()}
dists = {i:0 for i in [*positions]}
for t in range(100):
    mod.walk(mode = 'main', T = 1, verbose = True)
    new_positions = {i:(v.x, v.y) for i, v in mod.comm[0].items()}
    print('########')
    print('tstep= %i' % t)
    for i,v in new_positions.items():
        dist0 = 0
        distnot0 = 0
        if i in [*positions]:
            dist =np.sqrt((positions[i][0] - new_positions[i][0])**2 + (
                positions[i][1] - new_positions[i][1])**2)
            dists[i] = dist
            if dist != 0:
                print('ind  = %i1' % i)
                print('dist = %0.2f' % dist)
                print('old p= (%0.2f, %0.2f)' % (positions[i][0], positions[i][1]))
                print('new p= (%0.2f, %0.2f)' % (new_positions[i][0], new_positions[i][1]))
                distnot0+=1
            else:
                dist0 +=1
    print('ndist0   = %i\nndistnot0=%i\n~~~~~~~~~' % (dist0, distnot0))
    positions = new_positions
