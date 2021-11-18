library(shiny)
library(GWASmeta)

ui <- fluidPage(
  titlePanel("SMetABF for GWAS Meta-analysis"),
  fluidRow(
    column(4,
      p("Please upload the corresponding format of file. Example data format:"),
      tabsetPanel(
        tabPanel("single SNP",img(src = "single.png")),
        tabPanel("multiple SNPs",img(src = "multiple.png",height = 280, width = 450))
      ),
      #p(uiOutput("forma")),
      fileInput("upload","Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      selectInput("type", "How many SNPs do the data contain?", 
                  choices = c("single SNP","multiple SNPs")),
      checkboxInput("header", "Header", TRUE),
      selectInput("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      selectInput("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"')
      
    ),
    column(4,
        h4("parameters for calculation"),
        numericInput("num1", "prior sigma",
                       value = 0.5, min = 0, max = 1),
        radioButtons("prior", "Choose prior model", choices = c("fixed","correlated","independent")),
        numericInput("num2", "prior rho",
                       value = 0.5, min = 0, max = 1),
        sliderInput("iter", "iteration set",
                      value = 50, min = 0, max = 200),
        radioButtons("return", "Choose the return form", choices = c("log10","log2","origin")),
        actionButton("do", "Confirm")
    ),
    column(3,
      h2("results"),
      tags$hr(),
      tableOutput("contents")
    )
  )   
)

server = function(input,output,session){
  observeEvent(input$do,{
    output$contents <- renderTable({
      # input$upload will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      output$forma <- renderUI(input$type)
      prior.sigma <- input$num1
      prior.cor <- input$prior
      prior.rho <- input$num2  
      if(input$return=="log10"){
        log=FALSE
        log10=TRUE
      }else if(input$return=="log2"){
        log=TRUE
        log10=FALSE
      }else{
        log=FALSE
        log10=FALSE
      }
      n.iter <- input$iter
      req(input$upload)
      
      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
      tryCatch(
        {
          df <- read.csv(input$upload$datapath,
                         header = input$header,
                         sep = input$sep,
                         quote = input$quote)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
      if(input$type == "single SNP") {
        betas = df[,1]
        ses = df[,2]
        abf = shotgun.abf(betas,ses,prior.sigma,prior.cor,prior.rho,
                          cryptic.cor=NA,log,log10,na.rm,tolerance=1e-1000,n.iter,B=5)
      }
      else {
        get_abf <- function(i){
          SNP <- df[i,1]
          betas <- df[i,vbetas]
          ses <- df[i,vses]
          abfi <- shotgun.abf(betas,ses,prior.sigma,prior.cor,prior.rho,
                              cryptic.cor=NA,log,log10,na.rm,tolerance=1e-1000,n.iter,B=5)
          return(c(SNP,abfi))
        }
        vbetas <- seq(2,ncol(df),2)
        vses <- seq(3,ncol(df),2)
        re <- sapply(seq(1,nrow(df)),get_abf)
        abf <- data.frame(SNP=re[1,],ABF=re[2,])
      }
      return(abf)
    })
  })
}

shinyApp(ui, server)