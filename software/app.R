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
    column(3,
        h4("parameters for calculation"),
        numericInput("num1", "prior sigma",
                       value = 0.5, min = 0, max = 1),
        radioButtons("prior", "Choose prior model", choices = c("fixed","correlated","independent")),
        numericInput("num2", "prior rho",
                       value = 0.5, min = 0, max = 1),
        sliderInput("iter", "iteration set",
                      value = 50, min = 0, max = 500),
        radioButtons("return", "Choose the return form", choices = c("log10","log2","origin")),
        actionButton("do", "Confirm"),
        #p(uiOutput("forma")),
        p(uiOutput("forma1"))
    ),
    column(3,
      h2("top 20 lines of results"),
      downloadButton("downloadData", "Download"),
      tags$hr(),
      tableOutput("showData")
    )
  )   
)

server = function(input,output,session){
  observeEvent(input$do,{
    #output$forma <- renderText("The program is running, please wait...")
    calABF <- function(input,output){
      prior.sigma <- input$num1
      prior.cor <- input$prior
      prior.rho <- input$num2
      return <- input$return
      if(return=="log10"){
        log=FALSE
        log10=TRUE
      }else if(return=="log2"){
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
        ABF = shotgun.abf(betas,ses,prior.sigma,prior.cor,prior.rho,
                          cryptic.cor=NA,log,log10,na.rm,tolerance=1e-1000,n.iter,B=5)
        abf = data.frame(ABF=ABF)
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
        abf$ABF <- as.numeric(abf$ABF)
        abf$ABF <- round(abf$ABF,4)
        abf <- arrange(abf,desc(ABF))
      }
      return(abf)
    }
    abfL <- reactive({
      calABF(input,output)
    })
    abf <- abfL()
    output$forma1 <- renderText("Program has been finished")
    output$showData <- renderTable({
      head(abf,20)
    })
    output$downloadData <- downloadHandler(
      filename = function() {
        "abf_results.csv"
      },
      content = function(file) {
        write.csv(abf, file, row.names = FALSE)
      })
  })
}

shinyApp(ui, server)
