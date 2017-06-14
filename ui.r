library(shiny)
require(visNetwork)

shinyUI(fluidPage(
    titlePanel(title="Communication I4CS"),
    sidebarLayout(
      sidebarPanel(
        h3("Parameters:"),
        sliderInput("seed",
                    "Random generator seed:",
                    min = 1,
                    max = 100,
                    value = 50,
                    step = 1),
        sliderInput("dataselect",
                    "Observer:",
                    min = 1,
                    max = 50,
                    value = 25,
                    step = 1),
        sliderInput("x",
                    "x deviation:",
                    min = -2,
                    max = 2,
                    value = 0,
                    step = 0.1),
        sliderInput("y",
                    "y deviation:",
                    min = -2,
                    max = 2,
                    value = 0,
                    step = 0.1),
        sliderInput("alphaselect",
                    "Number of Neighbors:",
                    min = 2,
                    max = 50,
                    value = 10,
                    step = 1),
        sliderInput("top",
                    "Top of Special Cases:",
                    min = 2,
                    max = 50,
                    value = 3,
                    step = 1)
      ),
      mainPanel(
        h2("I4CS 2017 : Detection of special cases - Simulation"), 
        tabsetPanel(type="tab",
                    tabPanel("Read me", htmlOutput("readme") ),
                    tabPanel("Data Set", plotOutput("dataset") ),
                    tabPanel("Observations", plotOutput("ranking") ),
                    tabPanel("Neighborhoods", plotOutput("neighbors") ),
                    tabPanel("Neighbor", plotOutput("neighbor") ),
                    tabPanel("Rareness", plotOutput("rareness") ),
                    tabPanel("Special Case", plotOutput("cases") ),
                    tabPanel("Size of Neighborhoods", plotOutput("alpha") ),
                    tabPanel("Detection of Special Cases", plotOutput("result") )
        )
      )
    )
  )
)