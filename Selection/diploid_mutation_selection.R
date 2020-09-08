#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)

cwd <- 4

#user interface
ui <- fluidPage(pageWithSidebar( 
  
  headerPanel = headerPanel("Diploid model of selection with mutation"),
  
  sidebarPanel(
    
    HTML("<p style='font-size:14px'><B>Frequency of allele A over time (p).</B>"),

    selectInput(inputId = "u", label = "Mutation rate", choice = c(0,0.000001,0.00001,0.0001,0.001,0.01,0.1), 
                selected = 0.0),
  
    sliderInput(inputId = "wAA", label = HTML(paste0("W",tags$sub("AA"))), value = 1, 
                min = 0, max = 1, step = 0.01),
    sliderInput(inputId = "wAa", label = HTML(paste0("W",tags$sub("Aa"))), value = 0.9, 
                min = 0, max = 1, step = 0.01),
    sliderInput(inputId = "waa", label = HTML(paste0("W",tags$sub("aa"))), value = 0.81, 
                min = 0, max = 1, step = 0.01),
    
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
    plotOutput(outputId = 'viz'),
    
    #column(width = cwd,
    #checkboxInput(inputId = "w_plot",
    #              label = strong("Show average fitness by p"),
    #              value = FALSE)
    #),
    
    #column(width = cwd, 
    #checkboxInput(inputId = "delta_plot",
    #              label = strong("Show change in p by p"),
    #              value = FALSE)
    #),
    
    #column(width = cwd,
    #checkboxInput(inputId = "time_plot",
    #              label = strong("Show p by generation"),
    #              value = FALSE)
    #),
    
    selectInput("select", label = "Plot options", 
                choices = list("Average fitness by p" = 1, "Change in p by p" = 2,
                               "Allele frequency over time" = 3, "Genotypes over time" = 4), selected = 3)
    
    
  )
))

#back end code and response to user input
server <- function(input, output){
  
  output$viz <- renderPlot({
    
    p <- seq(0, 1, length.out = 1000)
    #parameters
    wAA = input$wAA
    wAa = input$wAa
    waa = input$waa
    p_0 = input$p_0
    u = as.numeric(input$u)
    gen = input$gen
    t = seq(1, gen, length.out = 1000)
    
    #if(input$w_plot){
    if(input$select == 1){
      W <- ((1-u)*p+u*(1-p))^2*wAA + 2*((1-u)*p+u*(1-p))*(u*p+(1-u)*(1-p))*wAa + (u*p+(1-u)*(1-p))^2*waa
      data.frame(p=p, W=W) %>%
        ggplot(aes(x = p, y = W)) +
        geom_line(color="firebrick", size=2) +
        xlab("Mean Fitness") +
        xlab("Allele frequency, p") +
        xlim(0, 1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
      
      #} else if(input$delta_plot){
    } else if(input$select == 2){
      W <- ((1-u)*p+u*(1-p))^2*wAA + 2*((1-u)*p+u*(1-p))*(u*p+(1-u)*(1-p))*wAa + (u*p+(1-u)*(1-p))^2*waa
      delta_p <- (((1-u)*p+u*(1-p))^2*wAA + ((1-u)*p+u*(1-p))*(u*p+(1-u)*(1-p))*wAa) / W - p
      data.frame(p=p, delta_p=delta_p) %>%
        ggplot(aes(x = p, y = delta_p)) +
        geom_line(color="firebrick", size=2) +
        geom_hline(yintercept = 0, lty = 2) +
        ylab(expression(paste(Delta,p))) +
        xlab("Allele frequency, p") +
        xlim(0, 1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
      
      #} else if(input$time_plot){
    } else if(input$select == 3){
      
      p_t <- rep(NA, gen)
      p_t[1] <- p_0
      for(i in 2:gen){
        W <- ((1-u)*p_t[i-1]+u*(1-p_t[i-1]))^2*wAA + 2*((1-u)*p_t[i-1]+u*(1-p_t[i-1]))*(u*p_t[i-1]+(1-u)*(1-p_t[i-1]))*wAa + (u*p_t[i-1]+(1-u)*(1-p_t[i-1]))^2*waa
        p_t[i] <- (((1-u)*p_t[i-1]+u*(1-p_t[i-1]))^2 * wAA + ((1-u)*p_t[i-1]+u*(1-p_t[i-1])) * (u*p_t[i-1]+(1-u)*(1-p_t[i-1])) * wAa) / W
      }
      
      data.frame(t = 1:gen, p_t = p_t) %>%
        ggplot(aes(x = t, y = p_t)) +
        geom_line(color="firebrick", size=2) +
        xlab("Generation") +
        ylab("Allele frequency, p") +
        ylim(0, 1) + 
        scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
    } else if(input$select == 4){
      
      p_t <- rep(NA, gen)
      p_t[1] <- p_0
      for(i in 2:gen){
        W <- ((1-u)*p_t[i-1]+u*(1-p_t[i-1]))^2*wAA + 2*((1-u)*p_t[i-1]+u*(1-p_t[i-1]))*(u*p_t[i-1]+(1-u)*(1-p_t[i-1]))*wAa + (u*p_t[i-1]+(1-u)*(1-p_t[i-1]))^2*waa
        p_t[i] <- (((1-u)*p_t[i-1]+u*(1-p_t[i-1]))^2 * wAA + ((1-u)*p_t[i-1]+u*(1-p_t[i-1])) * (u*p_t[i-1]+(1-u)*(1-p_t[i-1])) * wAa) / W
      }
      
      data.frame(t = 1:gen, p_t = p_t) %>%
        ggplot(aes(x = t)) +
        scale_color_manual(values = c("genotype AA" = "firebrick", "genotype Aa" = "forestgreen", "genotype aa" = "steelblue")) +
        geom_line(aes(y = (1-p_t)^2, color="genotype aa"), size=1) +
        geom_line(aes(y = 2*p_t*(1-p_t), color="genotype Aa"), size=1) +
        geom_line(aes(y = p_t^2, color="genotype AA"), size=1) +
        labs(x= "Generation",y="Allele frequency, p",
             color = "Legend") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) + 
        #        theme(panel.background =  
        #                element_rect(fill =  rgb(30, 144, 255, 25, 
        #                                         maxColorValue = 255)),
        #              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
    }
  })
}
# Run the application 
shinyApp(ui = ui, server = server)
