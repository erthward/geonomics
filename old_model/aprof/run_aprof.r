#First run all the code, to get a starting pop
source('../HART_final_project_no_rmd.r')


#Then dump the main model to a file
dump("core.model", "core_model.r"

#Then source that file
source("code_model.r")

#Create a tempfile
tmp = tempfile()

#Profile the core.model function
Rprof(tmp, line.profiling=T)
core.model(p.cross = 0.8, start.pops = START.POPS, n.time = N.TIME, size.pop = SIZE.POP, mean.d.move = MEAN.D.MOVE, for.avg = FOR.AVG, repro.cycle = REPRO.CYCLE, d.pair = D.PAIR, n.loci = N.LOCI, fit.advntg = FIT.ADVNTG, mean.off = MEAN.OFF, sd.off = SD.OFF, mean.disp.dist = MEAN.DISP.DIST, sd.disp.dist = SD.DISP.DIST)   
Rprof(append = F)

#Then create an aprof object
require(aprof)
a.o = aprof("core_model.r", tmp)

#Display info
a.o

#Summarize
summary(a.o)

#Plot (hard to see, but usable if I shrink the text of the core_model.r file and compare them; not much new information that can't be gleaned from the summary object)
plot(a.o)
