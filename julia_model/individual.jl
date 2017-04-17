module individual

export Individual, create_individual, functions


##########################################
#=

Module name:              individual


Module contains:
                          - definition of the Individual type
                          - function for creation of a simulated individual
                          - associated functions


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               12-28-15
Documentation:            URL


=#
##########################################



type Individual
  G::Genome #individual's x-ploid genotype
  x::Float64 #x coord
  y::Float64 #y coord

end



function sim_individual()

  #use genome.sim_genome and genome_arch variable to simulate individual's genome
  G = genome.sim_genome(genome_arch)A
  
  #randomly assign individual a valid starting location
  x =   #HOW TO SET?
  y =   #HOW TO SET?
  
  
  return Individual(G, x, y)

end

