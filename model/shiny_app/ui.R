library(shiny)

shinyUI(fluidPage(

    titlePanel('Landscape Genetics Model'),

    sidebarLayout(
        sidebarPanel(width = 2,
            helpText("Reset parameters to rerun model."),
            sliderInput("p.cross",
                        "Probability of crossing barrier:",
                        min = 0,
                        max = 1,
                        value = 0.5),
            numericInput("n.time",
                        "Number of timesteps:",
                        100),
            numericInput("n.loci",
                        "Number of loci to simulate:",
                         5),
            numericInput("r.dim",
                        "Raster size (cells per side):",
                         40),
            numericInput("size.pop",
                        "Size of starting populations:",
                         15),
            sliderInput("d.pair",
                        "Max distance between reproducers:",
                        min = 0.0005,
                        max = 0.5,
                        value = 0.005),
            selectInput("metric",
                        label = "Choose a metric to display",
                        choices = c('Gst', 'prop.h'),  #NOTE: CHANGE HERE TO INCLUDE VARIOUS METRICS ONCE READY!
                        selected = 'Gst')

            ),
    mainPanel(
        plotOutput("Fst.plot")
        )
       )
    ))


                        
