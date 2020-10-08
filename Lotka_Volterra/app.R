#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(phaseR)

lotkaVolterra <- function(t, y, parameters){
    x <- y[1]
    y <- y[2]
    lambda  <- parameters[1]
    epsilon <- parameters[2]
    eta     <- parameters[3]
    delta   <- parameters[4]
    dy    <- numeric(2)
    dy[1] <- lambda*x - epsilon*x*y
    dy[2] <- eta*x*y - delta*y
    list(dy)
}



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Modelos de Lotka-Volterra"),

    # Sidebar with a slider input for number of bins 
    tabsetPanel(
        tabPanel("Depredador-Presa",
                 sidebarLayout(
                     sidebarPanel(
                         numericInput("N1",
                                      "Número de presas (Nt):",
                                      min = 0,
                                      max = 100,
                                      value = 2),
                         numericInput("N2",
                                      "Número de depredadores (Pt):",
                                      min = 0,
                                      max = 60,
                                      value = 6),
                         sliderInput("Lambda",
                                     "Lambda",
                                     min = 0, 
                                     max = 5,
                                     value = 2,
                                     step = 0.1),
                         sliderInput("Epsilon",
                                     "Epsilon",
                                     min = 0, 
                                     max = 5,
                                     value = 1,
                                     step = 0.1),
                         sliderInput("Delta",
                                     "Delta",
                                     min = 0, 
                                     max = 5,
                                     value = 3,
                                     step = 0.1),
                         sliderInput("Gamma",
                                     "Gamma",
                                     min = 0, 
                                     max = 5,
                                     value = 2,
                                     step = 0.1),
                         checkboxInput("Nullclin",
                                       "Mostrar clinas nulas"),
                         sliderInput("Tiempos",
                                     "Número de tiempos a modelar:",
                                     min = 10, 
                                     max = 50,
                                     value = 10,
                                     step = 1),
                         numericInput("ylim",
                                      "Máximo poblacional",
                                      min = 0,
                                      max = 100,
                                      value = 10)),
                     
                     # Show a plot of the generated distribution
                     mainPanel(
                         uiOutput("ex3"),
                         uiOutput("ex4"),
                         plotOutput("distPlot"),
                         plotOutput("Field")
                     )
                 )),
        tabPanel("Competencia",
                 sidebarLayout(
                     sidebarPanel(
                         numericInput("Spp1",
                                      "N inicial especie 1",
                                      min = 0,
                                      value = 10,
                                      step = 5),
                         numericInput("r1",
                                      "Tasa de crecimiento de la especie 1",
                                      value = 0.9,
                                      step = 0.1),
                         numericInput("K1",
                                      "Capacidad de carga especie 1",
                                      min = 0,
                                      value = 100,
                                      step = 5),
                         numericInput("Alpha1",
                                      "Alpha de la especie 1 sobre 2",
                                      min = 0,
                                      value = 0.4,
                                      step = 0.1),
                         numericInput("Spp2",
                                      "N inicial especie 2",
                                      min = 0,
                                      value = 10,
                                      step = 5),
                         numericInput("r2",
                                      "Tasa de crecimiento de la especie 2",
                                      value = 0.9,
                                      step = 0.1),
                         numericInput("K2",
                                      "Capacidad de carga especie 1",
                                      min = 0,
                                      value = 140,
                                      step = 5),
                         numericInput("Alpha2",
                                      "Alpha de la especie 2 sobre 1",
                                      min = 0,
                                      value = 0.3,
                                      step = 0.1),
                         sliderInput("TiemposComp",
                                     "Número de tiempos a modelar:",
                                     min = 10, 
                                     max = 50,
                                     value = 10,
                                     step = 1),
                         numericInput("ylimComp",
                                      "Máximo poblacional",
                                      min = 0,
                                      value = 150),
                         checkboxInput("NullclinComp",
                                       "Mostrar clinas nulas")
                     ),
                     mainPanel(
                         uiOutput("Comp1"),
                         uiOutput("Comp2"),
                         plotOutput("Comp_Plot"),
                         plotOutput("FieldComp")
                     )
                 ))
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$Comp1 <- renderUI({
        withMathJax(
            helpText('Ecuación lotka-volterra Competencia 1
               $$\\dfrac{dx_1(t)}{dt} = r_1 x_1 - \\Big( 1 -\\Big( \\frac{x_1 + \\alpha_{1,2} x_2}{K_1}\\Big) \\Big) $$'))
    })
    output$Comp2 <- renderUI({
        withMathJax(
            helpText('Ecuación lotka-volterra Competencia 2
               $$\\dfrac{dx_2(t)}{dt} = r_2 x_2 - \\Big( 1 -\\Big( \\frac{x_2 + \\alpha_{2,1} x_2}{K_2}\\Big) \\Big) $$'))
    })
    
    output$ex3 <- renderUI({
        withMathJax(
            helpText('Ecuación lotka-volterra presa
               $$\\dfrac{dN(t)}{dt} = \\lambda N(t) - \\epsilon N(t)P(t)$$'))
    })
    
    output$ex4 <- renderUI({
        withMathJax(
            helpText('Ecuación lotka-volterra depredador
               $$\\dfrac{dP(t)}{dt} = \\delta N(t) P(t) - \\gamma P(t)$$'))
    })

    output$distPlot <- renderPlot({
        lotkaVolterra.numericalSolution <-  numericalSolution(lotkaVolterra, y0 = c(input$N1, input$N2), 
                                                              tlim = c(0,input$Tiempos), 
                                                              type = "one",
                                                              parameters = c(input$Lambda, input$Epsilon, input$Delta, input$Gamma), 
                                                              col = c("green", "orange"), 
                                                              ylim = c(0, input$ylim),
                                                              state =c("Presa", "Depredador"))
    })
    
    output$Comp_Plot <- renderPlot({
        comp.numericalSolution <-  numericalSolution(competition, y0 = c(input$Spp1, input$Spp2), 
                                                              tlim = c(0,input$TiemposComp), 
                                                              type = "one",
                                                              parameters = c(input$r1, input$K1, input$Alpha1, input$r2, input$K2, input$Alpha2), 
                                                              col = c("green", "orange"), 
                                                              ylim = c(0, input$ylimComp),
                                                              state =c("Especie 1", "Especie 2"))
    })
    
    output$Field <- renderPlot({
        lotkaVolterra.flowField <- flowField(lotkaVolterra, xlim = c(0, input$ylim), 
                                             ylim = c(0, input$ylim),
                                             parameters = c(input$Lambda, input$Epsilon, input$Delta, input$Gamma),
                                             points = 19, 
                                             add = FALSE,
                                             state.names = c("Presa", "Depredador"))
        grid()
        if(input$Nullclin){
            lotkaVolterra.nullclines <- nullclines(lotkaVolterra, 
                                                   xlim = c(-1, input$ylim), 
                                                   ylim = c(-1, input$ylim),
                                                   state.names = c("Presa", "Depredador"),
                                                   parameters = c(input$Lambda, input$Epsilon, input$Delta, input$Gamma), 
                                                   points = 20)   
        }
        
        y0 <- matrix(c(input$N1, input$N2), ncol = 2, nrow = 1, byrow = TRUE)
        lotkaVolterra.trajectory <- trajectory(lotkaVolterra, y0 = y0, tlim = c(0,input$Tiempos),parameters = c(input$Lambda, input$Epsilon, input$Delta, input$Gamma), col = rep("black", 1))
    })
    
    output$FieldComp <- renderPlot({
        
        High1 <- max(c(input$K1, input$K1/input$Alpha2))*1.05
        High2 <- max(c(input$K2, input$K2/input$Alpha1))*1.05
        
        Comp.flowField <- flowField(competition, xlim = c(0, High1),
                                    ylim = c(0, High2), 
                                             parameters = c(input$r1, input$K1, input$Alpha1, input$r2, input$K2, input$Alpha2),
                                             points = 50, 
                                             add = FALSE,
                                             state.names = c("Especie 1", "Especie 2"))
        grid()
        if(input$NullclinComp){
            lotkaVolterra.nullclines <- nullclines(competition, 
                                                   xlim = c(-1, High1), 
                                                   ylim = c(-1, High2),
                                                   state.names = c("Especie 1", "Especie 2"),
                                                   parameters = c(input$r1, input$K1, input$Alpha1, input$r2, input$K2, input$Alpha2), 
                                                   points = 19)   
        }
        
        #y0 <- matrix(c(input$N1, input$N2), ncol = 2, nrow = 1, byrow = TRUE)
        #lotkaVolterra.trajectory <- trajectory(lotkaVolterra, y0 = y0, tlim = c(0,input$Tiempos),parameters = c(input$Lambda, input$Epsilon, input$Delta, input$Gamma), col = rep("black", 1))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
