#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

library(tidyverse)

#Haploid
p_t <- function(t, p_0, WA, Wa, u){
  temp <- sqrt(4*u^2*Wa*WA+(1-u)^2*(WA-Wa)^2)
  temp2<-(((1-u)*(WA+Wa)-temp)^t*(p_0*temp-2*(1-p_0)*u*WA-p_0*(1-u)*(WA-Wa))+((1-u)*(WA+Wa)+temp)^t*(p_0*temp+2*(1-p_0)*u*WA+p_0*(1-u)*(WA-Wa)))
  temp3<-(((1-u)*(WA+Wa)-temp)^t*(temp+(1-2*p_0)*(1-2*u)*(WA-Wa)-u*(WA+Wa))+((1-u)*(WA+Wa)+temp)^t*(temp-(1-2*p_0)*(1-2*u)*(WA-Wa)+u*(Wa+WA)))
  temp2/temp3
}


#user interface
ui <- pageWithSidebar( 
  
  headerPanel = headerPanel("BIOL336: Haploid model of selection with mutation"),
  
  sidebarPanel(
    
    HTML("<p style='font-size:14px'><B>Frequency of allele A over time (p).</B>"),
    
    selectInput(inputId = "u", label = "Mutation rate", choice = c(0,0.000001,0.00001,0.0001,0.001,0.01,0.1), 
                selected = 0.0),
    
    sliderInput(inputId = "WA", label = HTML(paste0("W",tags$sub("A"))), value = 1, 
                min = 0, max = 2, step = 0.01),
    
    sliderInput(inputId = "Wa", label = HTML(paste0("W",tags$sub("a"))), value = 0.9, 
                min = 0, max = 2, step = 0.01),
    
    sliderInput(inputId = "p_0", label = "Initial p", value = 0.1, 
                min = 0, max = 1, step = 0.01),
    
    sliderInput(inputId = "gen", label = "Number of generations", value = 100, 
                min = 2, max = 200, step = 1),

    HTML("<p style='font-size:12px'>Take me to:
         <UL><LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/HaploidSelection/'>Haploid selection<a> (<A href='https://shiney.zoology.ubc.ca/otto/HaploidMutationSelection/'>with mutation<a>)
         <LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/DiploidSelection/'>Diploid selection<a> (<A href='https://shiney.zoology.ubc.ca/otto/DiploidMutationSelection/'>with mutation<a>)
         <LI style='font-size:12px'><A href='https://shiney.zoology.ubc.ca/otto/DiploidDriftSelection/'>Diploid selection with drift<a>
         </UL><P>"),
    
    HTML("<p style='font-size:8px'>Modified from: Copyright (c) 2017 Silas Tittes, MIT License, https://github.com/silastittes/shiny_popgen</p>")
    
  ), 
  
  mainPanel =  mainPanel(
    plotOutput(outputId = 'freq')
  )
)

#back end code and response to user input
server <- function(input, output){
  
  sim_A <- reactive({
    
    #parameters
    WA = input$WA
    Wa = input$Wa
    p_0 = input$p_0
    u = as.numeric(input$u)
    gen = input$gen
    t = seq(1, gen, length.out = 1000)
    
    return(data.frame(t=t, p = p_t(t, p_0, WA, Wa, u)))
    
  })
  
  
  output$freq <- renderPlot({
    
    sim_A() %>%
      ggplot(aes(x = t, y = p)) +
      geom_line(color="firebrick", size=2) +
      xlab("Generations") +
      ylab("Allele frequency, p") +
      ylim(0, 1)+ 
      scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) + 
      # theme(panel.background =  
      #           element_rect(fill =  rgb(30, 144, 255, 25, 
      #                                    maxColorValue = 255)),
      #       text = element_text(size=16, family= "Times"))+ 
      theme_classic()+theme(text = element_text(size=20))
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
