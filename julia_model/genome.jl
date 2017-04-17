module genome

export Chromosome, Genome, rand_bern, sim_G, set_l_c, define_dist_s, define_dist_D, sim_s, sim_D, build_genomic_arch, sim_genome 

using Distributions


##########################################
#=

Module name:              genome

Module contents:          - definition of the Genome type
                          - function for generation of a simulated genome
                          - associated functions


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               12-28-15
Documentation:            URL


=#
##########################################





#----------------------------------
# TYPES ---------------------------
#----------------------------------


#type Chromosome
#  l::Int # length (i.e. number of markers)
#  G::Array{Int} # genotype (i.e. array of alleles, biallelic, i.e. binary)
#  s::Array{Float16} #array of selection coefficients for all loci
#  D::Array{Float16} #array of map distances (i.e. physical linkage disequilibria) between each locus and the following
#  #NOTE: This implies that the last values in this array must be 0!
#end
     

type Genomic_architecture
  x::Int #ploidy (NOTE: for now will be 2 by default; later could consider enabling polyploidy)
  n::Int #haploid number of chromosomes
  L::Int #total length (i.e. number of markers)
  l_c::Array{Int} #length of each chromosome
  p::Dict #Dict of dominant (i.e. coded by '1', not '0') allele frequencies for all (numbered by keys) chromosomes in haploid genome
  s::Dict #Dict of selection coefficients, for all chroms
  D::Dict #Dict of linkage disequilibria (i.e. map distances) between each locus and the next, for all chroms (NOTE: because of this definition, the last value is forced to 0)
  sex::Bool
end





type Genome
  G::Dict #Dict of numbered i*j arrays, each subarray containing the biallelic genotype data for all j (for a diploid, 2) copies of chromosome i

  #NOTE: should/could I create methods herein for mutation, crossing over, etc??

#NOTE: need to ensure that size(G)[1] == n and that size of each subarray in G == x*l_c

end









#----------------------------------
# FUNCTIONS -----------------------
#----------------------------------


#function for choosing a random number from a Bernoulli distribution of an arbitrary p value
#(to be mapped to allele frequencies in order to generate genotypes)
function rand_bern(prob::Float16)
  return rand(Bernoulli(prob))
end



#generate allele_freqs
function gen_allele_freqs(l::Int)
  return rand(Beta(1,1),l)
end



#simulate genotypes, based on an array of major allele freqs
function sim_G(p::Array) 
  return map(rand_bern, p)
end


#define distribution for selection coefficients
function define_dist_s(args...; alpha_s = 0.15, beta_s = 2)
#function define_dist_s(alpha_s::Float64, beta_s::Float64)
  if false #NOTE: replace here with "if certain command-line flag provided, define distribution manually according to arguments provided
    return Union{} #replace this with ability to call whatever distribution requested, using args... from above
  else
    return Beta(alpha_s, beta_s) 
  end
  #NOTE: For now, beta seems an intuitive and fine way to model this
  #For effectively no selection: Beta(1,1e5)
  #For spread of different selection values, but highly zero-inflated: ~ Beta(0.15, 2)
  #For a wide range of selection coefficients, something like: Beta(1,5)
  #For a uniform range of selection coefficients between 0 and 1: beta (1,1)
end


#simulate selection coefficients
function sim_s(l::Int, dist_s::Distributions.Distribution)
  return rand(dist_s, l)
end

#define distribution for physical linkage disequilibria (i.e. map distances)
function define_dist_D(args...; alpha_D = 3e3, beta_D = 7e3)
  if false #NOTE: replace here with "if certain command-line flag provided, define distribution manually according to arguments provided
    return Union{} #replace this with ability to call whatever distribution requested, usings args... from above
  else
    return Beta(alpha_D, beta_D)
    #NOTE: for now, using the Beta, which can be very flexibly parameterized
    #NOTE: current default alpha/beta vals, after being subtracted from 0.5 in sim_D function, will result in a tight distribution of D vals around 0.21 (~ 95% between 0.19 and 0.22)
    #NOTE: to essentially fix D at 0.5, Beta(1,1e7) will do...
    #NOTE: Ideally, would be good to look into to developing some sort of mixture distribution to reasonably and flexibly model map-distance values...
  end
end


#simulate linkage values
function sim_D(l::Int, dist_D::Distributions.Distribution)
  return 0.5 .- rand(dist_D, l) 
end


#set chromosome lengths
function set_l_c(L, n, args... ; even_chrom_sizes = true)
  if !even_chrom_sizes
    l_c = Union{} #NOTE: Instead, if command-line arg provided, allow it to use user-supplied array stipulating of chromosome lengths
    return l_c
  else
    l_c = [Int(round(L/n)) for i in 1:n]
    return l_c
  end
end



#build the genomic architecture
#NOTE: This will create the "template" for the genomic architecture of the hapolid genome that will then be used to simulate individuals and populations
function build_genomic_arch(L::Int, n::Int ; x = 2, sex = false) 
         #NOTE: x = ploidy, for now set to 2 (i.e.  diploidy)
         #NOTE: how to operationalize sexuality?! for now defaults to false
  l_c = set_l_c(L, n)
  dist_s = define_dist_s()
  dist_D = define_dist_D()
  p = Dict()
  s = Dict()
  D = Dict()
  for chrom in 1:n
    p[chrom] = gen_allele_freqs(l_c[chrom])
    s[chrom] = sim_s(l_c[chrom], dist_s)
    D[chrom] = sim_D(l_c[chrom], dist_D)
    D[chrom][end] = 0 #because each D value expresses linkage between that locus and the next, last value is forced to 0!
  end
  return Genomic_architecture(x, n, sum(l_c), l_c, p, s, D, false), dist_s, dist_D
end




#simulate genome
function sim_genome(genomic_arch)
  genome = Dict()
  for chrom in 1:genomic_arch.n
    chromosome = ones(genomic_arch.l_c[chrom], genomic_arch.x)*999 #if for some reason any loci are not properly set to either 0 or 1, they will stick out as 999's
    for copy in 1:genome_arch.x
      chromosome[:,copy] = sim_G(genomic_arch.p[chrom])
    end

    genome[chrom] = chromosome
  end
  return Genome(genome)
end




end #of module
