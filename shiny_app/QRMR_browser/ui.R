library(shiny)
library(plotly)
library(ggplot2)
library(shinythemes)
`%>%` <- magrittr::`%>%`


# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    theme = shinytheme("spacelab"),
    # Application title
    titlePanel("QRMR Browser"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            h2('About'),
            p('This app displays output from quantile regression + meta regression on UK Biobank genotyping data. Three quantitative traits were tested: body mass index, height, and hemoglobin A1c. Quantile regression and standard errors were determind for quantiles 0.1, 0.25, 0.5, 0.75, and 0.9. Meta regression was then performed on these quantile estimates using a mixed-effects model represented by the equation below:'),
            
            uiOutput("formula"),
            span('For more information on methods, check out the '), 
            a('slide deck', href = 'https://brandonsie.github.io/pages/QR_slides.html'),
            span('.'),
            
            h2('Data Filters'),
            column(6, 
            sliderInput("b0p_range",
                        "log10p_b0:",
                        min = 0,
                        max = 306,
                        value = c(0,306)),
            sliderInput("b1p_range",
                        "log10p_b1:",
                        min = 0,
                        max = 307,
                        value = c(0,307)),
            sliderInput("b2p_range",
                        "log10p_b2:",
                        min = 0,
                        max = 303,
                        value = c(0,303))
            ), 
            column(6,
            sliderInput("b0_range",
                        "b0_estimate:",
                        min = -17,
                        max = 8,
                        value = c(-17,8)),
            sliderInput("b1_range",
                        "b1 estimate:",
                        min = -18,
                        max = 26,
                        value = c(-18,26)),
            sliderInput("b2_range",
                        "b2 estimate:",
                        min = -80,
                        max = 73,
                        value = c(-80,73))
            ),
            selectInput("phenotype", "Choose Phenotype:", as.list(c("bmi", "ht", "a1c"))),
            uiOutput("choose_rsid"),
            uiOutput("choose_phenotype"),
            
            h2('Figures'),
            p('QRMR coefficient scatterplot provides a scatterplot of b1 and b2 coefficients from meta-regression, faceted by phenotype, subset to variants with significant non-homogenous effect size after multiple hypothesis correction. Point color indicates whether or not a point matches the parameters selected in the Data Filters input portion of the left sidebar (above).'),
            p("Meta Regression effect size plot shows data for one SNP, with mean linear regression in red, quantile regression in black, and meta-regression in blue ribbon."),
            p("Simlulated Quantile-Dependent Distribution Shift shows the effect per copy of risk allele of the estimated quantile regression fit on a normal distribution with the real mean and standard deviation of that phenotype from the observed samples. The normal distribution may not accurately represent the real data distribution and is intended to provide a reference to help interpret patterns from quantile regression effect sizes."),
            p("The bottom table shows meta regression information for variants that match the parameters selected in the Data Filters input portion of the left sidebar (above).")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotlyOutput("beta_scatter"),
            # plotlyOutput("g_manhattan"),
            # tableOutput("plotdf"),
            plotOutput("g_qrmr"),
            tableOutput("table"),
            textOutput("this_qrmr")
        )
    )
))
