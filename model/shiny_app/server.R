#load necessary pacakges
library(shiny)
library(raster)
library(rgeos)
library(kriging)
library(spatstat)
library(circular)
library(poppr)
library(mmod)
library(adegenet)
#library(Geneland)
#library(ggplot2)

source('model_stuff.r')



#Run in 'debugging mode'?
#debug = T
debug = F



#Create running plot (T) or output PDF of plot (F)?
running.plot = T
#running.plot = F

#Create PDF (T) or plot to active device (F)?  (Only utilized if running.plot = F)
#PDF = T
PDF = F


########
#Colors for points
invisible(palette(c('#f03b20', '#ffeda0', '#feb24c'))); invisible(dev.off())
#palette(c('#2710DF', '#F6FF3D', '#60FFAE')); dev.off() 


#Colors for raster
Blues = c('#deebf7', '#9ecae1', '#3182bd')
ramp = colorRampPalette(Blues, space = 'Lab', bias = 0.5)


##############
#Plotting params
PCH = 21
FG = '#4C4C4C'
CEX = 1.2


#num.loci = c(5,10,15,20, 50)





#R.DIM = 40
#SIZE.POP = 15
#N.LOCI = 10 
#N.TIME = 100
#D.PAIR = 0.005
REPRO.CYCLE = 5
FOR.AVG = 6 
MEAN.OFF = 0.5
SD.OFF = 0.1
FIT.ADVNTG = 1.0
MEAN.DISP.DIST = 7
MEAN.D.MOVE = 0.01





###################################################################

####SHINY SERVER!
shinyServer(
   function(input, output) {


   output$Fst.plot <- renderPlot({


    p.cross = c(input$p.cross)

    N.TIME = input$n.time
    N.LOCI = input$n.loci
    R.DIM = input$r.dim
    SIZE.POP = input$size.pop
    D.PAIR = input$d.pair

    metric = input$metric



    set.seed(105)



    #create landscape layer 1 (heterogeneously distributed resource layer)
    pts <- data.frame(cbind(xCoord = runif(25,0,R.DIM), yCoord = runif(25,0,R.DIM), habitat = runif(25, 0, 1)))
    env <- kriging(x=round(pts$xCoord), y=round(pts$yCoord), response = pts$habitat)
    env <- 0.7*rasterFromXYZ(env$map)
    env@data@values = env@data@values/max(env@data@values) #normalize to 0 < cell_value < 1
    invisible(setExtent(env, extent(0,100,0,100)))
    
    
    
    
    #calculate mean.d.move from env and MEAN.D.MOVE
    MEAN.D.MOVE = MEAN.D.MOVE*sqrt((env@extent@xmax-env@extent@xmin)^2 + (env@extent@ymax-env@extent@ymin)^2)
    D.PAIR = D.PAIR*sqrt((env@extent@xmax-env@extent@xmin)^2 + (env@extent@ymax-env@extent@ymin)^2)
    
    
    
    
    hab.x = sort(runif(24,env@extent@xmin, env@extent@xmax))
    hab.y = sort(runif(24,env@extent@ymin, env@extent@ymax))
    hab.x = c(env@extent@xmin, env@extent@xmin, hab.x, env@extent@xmax, env@extent@xmax, env@extent@xmin)
    hab.y = c(env@extent@ymin, hab.y[1], hab.y, hab.y[length(hab.y)], env@extent@ymin, env@extent@ymin)
    hab.poly = SpatialPolygons(list(Polygons(list(Polygon(matrix(cbind(hab.x, hab.y), ncol = 2))), c('polygon'))))
    hab = rasterize(hab.poly, env)  #NOTE: let 1's = hab A (i.e. pop 'p1') and 0's = B (i.e. pop 'p2')
    hab@data@values[which(is.na(hab@data@values), arr.ind = TRUE)] = 0






    #Create timestep vector (for plotting output)
    TIME = c(1,seq(1:N.TIME)[which(seq(1:N.TIME) %% 5 ==0)])





    START.POPS = sim.pop(n.loci = N.LOCI, size.pop = SIZE.POP, for.avg = FOR.AVG, hab = hab, env = env)


    START.POPS = burn.in(start.pops = START.POPS, size.pop = SIZE.POP, mean.d.move = MEAN.D.MOVE, for.avg = FOR.AVG, hab = hab, env = env)




#results = list()



#for (p in p.cross){

  results =                               core.model(p.cross = p.cross,
  #results[[as.character(p.cross)]] = core.model(p.cross = p.cross,
                                          start.pops = START.POPS, 
                                          n.time = N.TIME, 
                                          size.pop = SIZE.POP, 
                                          mean.d.move = MEAN.D.MOVE, 
                                          for.avg = FOR.AVG, 
                                          repro.cycle = REPRO.CYCLE, 
                                          d.pair = D.PAIR,
                                          n.loci = N.LOCI, 
                                          fit.advntg = FIT.ADVNTG,
                                          mean.off = MEAN.OFF,
                                          sd.off = SD.OFF,
                                          mean.disp.dist = MEAN.DISP.DIST, 
                                          sd.disp.dist = SD.DISP.DIST, 
                                          hab = hab,
                                          env = env)

#}


#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

par(mfrow = c(1,3))

plot(env, col = ramp(5), main = list('Starting locations'), cex.lab = CEX, cex.axis = CEX)
lines(hab.x, hab.y, lwd = 1.5)
#points(results[[as.character(p.cross[1])]][['start']][, c('x', 'y')], bg = results[['start']]$pop, pch = PCH, fg = FG, cex = CEX)
start.pops = rbind(START.POPS[['p1']][[1]], START.POPS[['p2']][[1]])
points(start.pops[, c('x', 'y')], bg = start.pops$pop, pch = PCH, fg = FG, cex = CEX)


plot(env, col = ramp(5), main = list('Final locations'), cex.lab = CEX, cex.axis = CEX)
lines(hab.x, hab.y, lwd = 1.5)
#points(results[[as.character(p.cross[1])]][['end']][, c('x', 'y')], bg = results[['end']]$pop, pch = PCH, fg = FG, cex = CEX)
points(results[['all.pops']][, c('x', 'y')], bg = results[['all.pops']]$pop, pch = PCH, fg = FG, cex = CEX)


#plot(TIME, results[[as.character(p.cross[1])]][['metrics']][[metric]], type = 'l', col = 'red', lwd = 1.25,
plot(TIME, results[['metrics']][[metric]], type = 'l', col = 'red', lwd = 1.25,
           xlab = 'model time (timesteps)', ylab = metric, cex.lab = CEX, cex.axis = CEX, ylim = c(0,1),
           main = list(sprintf('%s as a function of model time', metric), cex = CEX))
#legend('bottomright', legend = as.character(p.cross), 
#         #lty = seq(length(p.cross)), 
#         lwd = rep(1.25, length(p.cross)),  
#         col = 'black', cex = CEX)




### CLOSE SHINY SERVER STATEMENT!
  })
})
