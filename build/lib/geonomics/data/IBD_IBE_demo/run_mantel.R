library(vegan)
library(phytools)

# load MMRR
source('./MMRR.R')

# function to read a headerless CSV of numerics in cwd with name <file_basename>.csv
read_matrix = function(file_basename){
    mat = as.matrix(read.table(paste0('./', file_basename, '.csv'), sep=','))
    return(mat)
}

# read the gen, env, and geo data
gen = read_matrix('gen')
env = read_matrix('env')
geo = read_matrix('geo')

# run the partial Mantel test and print the result
mod = mantel.partial(gen, env, geo)
print('PARTIAL MANTEL')
print(mod)

# run the opposite partial Mantel and print the result
mod = mantel.partial(gen, geo, env)
print('INVERSE PARTIAL MANTEL')
print(mod)

# run a multi-Mantel and print the result
mod = multi.mantel(gen, list(env, geo))
print('MULTI-MANTEL')
print(mod)

# run an MMRR model (results should be identical to the multi-Mantel)
mod = MMRR(gen, list(env, geo))
print('MMRR')
print(mod)
