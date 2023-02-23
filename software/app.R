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
      selectInput("type", "Single or multiple SNPs to be analyzed?",
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
    ),
    column(4,
      h2("top 20 lines of results"),
      downloadButton("downloadData", "Download"),
      tags$hr(),
      tableOutput("showData")
    )
  ),
  tags$head(
    tags$style(
      HTML(".shiny-notification {
             position:fixed;
             top: calc(62%);;
             left: calc(35%);;
             }"
      )
    )
  )
)

server = function(input,output,session){
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
      ABFl = shotgun_abf_model(betas,ses,prior.sigma,prior.cor,prior.rho,
                              cryptic.cor=NA,log,log10,na.rm,tolerance=1e-1000,n.iter,B=5)
      ABF = ABFl$ABF
      submodel = ABFl$model
      abf = data.frame(ABF=ABF,model=submodel)
    }
    else {
      vname <- 1
      vbetas <- seq(2,ncol(df),2)
      vses <- seq(3,ncol(df),2)
      nstud <- length(vbetas)
      if(length(vbetas)!=length(vses)){
        return("Vectors betas and ses do not have the same length!")
      }
      if(length(vbetas)==1){
        return("Only one study involved!")
      }
      get_abf <- function(i){
        SNP <- df[i,vname]
        betas <- df[i,vbetas]
        ses <- df[i,vses]
        nstudies <- df[i,"counts"]
        studiesUsed <- paste(1-as.numeric(is.na(betas)),collapse="")
        tryCatch(
          {
            abfL <- shotgun_abf_model(betas,ses,prior.sigma,prior.cor,prior.rho,
                                     cryptic.cor=NA,log,log10,na.rm,tolerance=1e-1000,n.iter,B=5)
            abfvalue <- abfL$ABF
            submodel <- abfL$model
          },
          error = function(e){
            abfvalue <- NA
            submodel <- "NA"
          })
        return(c(SNP,abfvalue,submodel,nstudies,studiesUsed))
      }
      df[df==0] <- NA
      df[df==Inf] <- NA
      for(i in 1:nstud){
        df[which(is.na(df[vbetas[i]])),vses[i]] <- NA
        df[which(is.na(df[vses[i]])),vbetas[i]] <- NA
      }
      df$counts <- nstud-rowSums(is.na(df[,vbetas]))
      re <- sapply(seq(1,nrow(df)),get_abf)
      abf <- data.frame(SNP=re[1,],ABF=re[2,],model=re[3,],n_studies=re[4,],studies_involved=re[5,])
      abf$ABF <- as.numeric(abf$ABF)
      abf$ABF <- round(abf$ABF,4)
      abf <- arrange(abf,desc(ABF))
    }
    return(abf)
  }
  observeEvent(input$do,{
    progress <- Progress$new(session, min=0, max=10)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    abfL <- reactive({
      calABF(input,output)
    })
    abf <- abfL()
    if(class(abf)=="data.frame"){
      output$showData <- renderTable({
        head(abf,20)
      })
    }else{
      output$showData <- renderTable({
        abf
      })
    }
    output$downloadData <- downloadHandler(
      filename = function() {
        "abf_results.csv"
      },
      content = function(file) {
        write.csv(abf, file, row.names = FALSE)
      })
    # withProgress(message = "Calculation in progress",
    #              detail = "This may take a while...", value = 0, {
    #                for (i in 1:30) {
    #                  incProgress(1/30)
    #                  Sys.sleep(0.25)
    #                }
    #             })
  })
}

shinyApp(ui, server)
shiny::runApp()