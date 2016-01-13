module population

export Population, create_population, functions

##########################################
#=

Module name:              population


Module contains:
                          - definition of the Population type
                          - function for creating a population of Individuals (incl. their genomes and associated data)
                          - associated functions


Author:                    Drew Ellison Hart
Email:                     drew.hart@berkeley.edu
Github:                    URL
Start date:                12-28-15
Documentation:             URL


=#
##########################################


type Population
  inds::Dict #Dict of all individuals' locations and genotypes
  birth_rate::Float64
  death_rate::Float64
                 
 
  #Add other demographic parameters??


function sim_population(N)
  pop = Dict()
  for i in range(1,N)
  # use individual.sim_individual to simulate individuals and add them to the population
    pop[i] = individual.sim_individual()

  return Population(pop, birth_rate, death_rate)

end

