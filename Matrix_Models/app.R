#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(igraph)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Modelo matricial"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("P1",
                        "Probabilidad de sobrevivir de huevo a juvenil:",
                        min = 0,
                        max = 1,
                        value = 0.2,
                        step = 0.05),
            sliderInput("P2",
                         "Probabilidad de sobrevivir de juvenil a adulto:",
                         min = 0,
                         max = 1,
                         value = 0.5,
                         step = 0.05),
            numericInput("F3",
                         "Fertilidad de Adultos:",
                         min = 0,
                         max = 1000,
                         value = 9),
            numericInput("N1",
                         "Población huevos (Etapa 1):",
                         min = 0,
                         max = 1000,
                         value = 20,
                         step = 10),
            numericInput("N2",
                         "Población jueveniles (Etapa 2):",
                         min = 0,
                         max = 1000,
                         value = 20,
                         step = 10),
            numericInput("N3",
                         "Población adultos (Etapa 3):",
                         min = 0,
                         max = 1000,
                         value = 20,
                         step = 10),
            numericInput("Years",
                         "Años a proyectar:",
                         min = 0,
                         max = 1000,
                         value = 50,
                         step = 5)
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Modelo",
                         plotOutput("distPlot"),
                         fluidRow(
                             column(2,tableOutput("Matrix")),
                             column(6,tableOutput("Ns"), offset = 1)
                         ),
                         uiOutput("Math")
                ),
                tabPanel("Projección",
                         plotOutput("Popplot"),
                         plotOutput("Elasticity")),
                tabPanel("Parámetros",
                         textOutput("lambda"),
                         plotOutput("eee"),
                         plotOutput("Repro"))
                
            )
        )
        
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    DF <- shiny::reactive({
        DF <- data.frame(from = c(1,2,3,3), 
                         to = c(2,3,3,1),  
                         S = c(input$P1, input$P2, 0, 0),
                         Fe = c(0,0,0,input$F3)) %>% 
            mutate(Scaled_S = S/max(S),
                   Scaled_Fe = Fe/max(Fe))
        
        Graph <- graph_from_data_frame(DF)
        list(DF = DF, Graph = Graph)
    })
    
    PopIn <- shiny::reactive({
        c(input$N1, input$N2, input$N3)
    })
    
    Matrix <- shiny::reactive({
        Mat <- matrix(data = 0, ncol = 3, nrow = 3)
        Mat[1,3] <- DF()$DF$Fe[4]
        Mat[2,1] <- DF()$DF$S[1]
        Mat[3,2] <- DF()$DF$S[2]
        colnames(Mat) <- c("", "","")
        Mat
    })
    
    output$Ns <- renderTable({
        data.frame(N = c(input$N1, input$N2, input$N3))
    }, colnames = F)

    output$distPlot <- renderPlot({
        ggraph(DF()$Graph, layout = 'linear')+
            geom_edge_link(aes(alpha = Scaled_S), color = "blue", arrow = arrow()) +
            annotate(x = 1.5 , y= 0.052, geom = "text", label= DF()$DF$S[1])+
            annotate(x = 2.5 , y= 0.052, geom = "text", label= DF()$DF$S[2])+
            annotate(x = 2 , y= -0.7, geom = "text", label=DF()$DF$Fe[4])+
            geom_edge_arc(aes(alpha = Scaled_Fe), color = "red", arrow = arrow()) +
            geom_node_point(size = 4)  +
            geom_node_label(aes(label = name)) +
            ggraph::theme_graph() + ggraph::scale_edge_alpha_continuous(range =c(0,1), guide = "none")
    })
    output$Matrix <- renderTable({
        Mat <- matrix(data = 0, ncol = 3, nrow = 3)
        Mat[1,3] <- DF()$DF$Fe[4]
        Mat[2,1] <- DF()$DF$S[1]
        Mat[3,2] <- DF()$DF$S[2]
        colnames(Mat) <- c("", "","")
        Mat
    }, colnames = F)
    
    output$Popplot <- renderPlot({
        result.matrix <-matrix(ncol = 3)
        result.matrix[1,]<- PopIn()
        
         for (t in 2:input$Years){
             pop.year<- result.matrix[t-1,]%*% Matrix()
             result.matrix<-rbind(result.matrix, pop.year)
             message(t)
        }
         tot.pop<-rowSums(result.matrix)
         matrix_boba<-cbind(result.matrix,tot.pop)
         
         colnames(matrix_boba)<-c("Huevo", "Juvenil", "Adulto", "Total")
         rownames(matrix_boba)<-seq(1,input$Years,1)
         df_bobo <- matrix_boba %>% as.data.frame() %>% 
             mutate(t = dplyr::row_number()) %>% 
             pivot_longer(-t, names_to = "Etapa", values_to = "n")
        ggplot(df_bobo, aes(x = t, y = n)) + geom_path(aes(color = Etapa)) + theme_bw()
    })
    
    output$lambda <- renderText({
        proj.matrix <- Matrix()
        # Paso 1. Calcula los eigenvalues y eigenvectors para la matriz de proyección
        eigs <- eigen(proj.matrix)
        
        dom.index <- which.max(eigs[["values"]])
        lambda <- Re(eigs[["values"]][dom.index])
        paste("Valor de","lambda", "=" ,round(lambda, 2))
    })
    
    output$eee <- renderPlot({
        proj.matrix <- Matrix()
        # Paso 1. Calcula los eigenvalues y eigenvectors para la matriz de proyección
        eigs <- eigen(proj.matrix)
        dom.index <- which.max(eigs[["values"]])
        w <- Re(eigs[["vectors"]][, dom.index])
        eee <- w / sum(w)
        
        # Step 2. Plot the stable stage distribution
        barplot(eee, las=1, names=c("Huevo", "Juvenil", "Adulto"), 
                main="Distribución Estable de Edades",
                xlab="Clase de edades", ylab="Proporción relativa")
        
    })
    
    output$Repro <- renderPlot({
        m <- eigen(t(Matrix()))
        v <- Re(m$vectors[, which.max(Re(m$values))])
        rv <- v/v[1]
        
        # Paso 2. Graficar el valor reproductivo
        barplot(rv, las=1, names=c("Huevo", "Juvenil", "Adulto"),
                xlab="Clase de edad", ylab="Valor reproductivo", main="Valor reproductivo por clase")
    })
    
    output$Elasticity <- renderPlot({
        
        proj.matrix <- Matrix()
        
        eigs <- eigen(proj.matrix)
        
        # Paso 2. Calcula la tasa finita de crecimiento (λ)
        dom.index <- which.max(eigs[["values"]])
        lambda <- Re(eigs[["values"]][dom.index])
        
        
        
        # Paso 1. Calcular EEE
        w <- Re(eigs[["vectors"]][, dom.index])
        eee <- w / sum(w)
        
        # Paso 1. Calcular el valor reproductivo
        m <- eigen(t(proj.matrix))
        v <- Re(m$vectors[, which.max(Re(m$values))])
        rv <- v/v[1]
        
                # Paso 1. Calcular sensitividad
        # (v y w los calculamos anteriormente)
        vw <- v %*% t(w)
        sens <- vw/as.numeric(v %*% w)
        
        # Paso 2. Calcular elasticidad
        # (lambda y sensitividad fueron calculadas anteriormente)
        elas <- (proj.matrix/lambda) * sens
        
        DF <- data.frame(Etapa = c("Huevo", "Juvenil", "Adulto"), Fertilidad = elas[1,], Supervivencia = c(unlist(sapply(1:2, function(x) elas[x+1, x])),0)) %>% 
            pivot_longer(cols =-Etapa, names_to = "Parametro", values_to = "Elasticidad")
        
        ggplot(DF, aes(x = Etapa, y = Elasticidad)) + geom_col(aes(fill = Parametro), position = "dodge") +
            theme_bw()
    }) 
}

# Run the application 
shinyApp(ui = ui, server = server)
