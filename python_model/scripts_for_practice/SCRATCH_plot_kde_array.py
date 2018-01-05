def plot_movement_surf_vectors(land):
    for i in range(len(land.movement_surf)):
        for j in range(len(land.movement_surf[0])):
            kde = land.movement_surf[i][j]
            samps = kde.sample(15)
            avg_dirs = []
            for samp in samps:
                avg_x = mean([np.cos(d)*sqrt(2) for d in samp])
                avg_y = mean([np.sin(d)*sqrt(2) for d in samp])
                unit_vector_avg = np.arctan(avg_y/avg_x)
                avg_dirs.append(unit_vector_avg)
    
            avg_avg = mean(avg_dirs)
    
            x = j  
            y = i
            #NOTE: would add 0.5 to the point above, to plot arrow in center of cell, but then would subtract 0.5 to account for visual reconciliation of offset between axes and raster, so net result is to use the cell i,j indices as the plotting point to root each arrow in its cell center
            dx = np.cos(avg_avg) * np.sqrt(2) / 2
            dy = np.sin(avg_avg) * np.sqrt(2) / 2
            plt.arrow(x, y, dx, dy, alpha = 0.75, color = 'black', head_width = 0.25, head_length = 0.25)


