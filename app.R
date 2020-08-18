#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Michael Sierk
# 7/29/20
# For Nephros, Inc.
# 
# Read in .csv from EPI2ME software, count reads from each taxid. Output:
#   1) Table of counts
#   2) Histograms

# csv headers
# filename	read_id	exit_status	runid	taxid	barcode	accuracy	lca	genus_taxid
# Allow user to select file from filesystem
library(shiny)
library(dplyr)
library(taxonomizr) # converts taxid to genus & species
library(ggplot2)
library(scales) # apparently needed for labels=comma below?
library(lessR) # for writing Excel output
library(readxl)
library(writexl)
source("getDBs.R") # script for updating DBs

# Change to TRUE to (re)generate the the taxonomy database nameNode.sqlite
make_taxonomy_db = FALSE
if (make_taxonomy_db) {
    makeTaxonomyDB()
}


# Define UI for application 
ui <- fluidPage(

    # Application title
    titlePanel("Genus Counts"),

    # Sidebar with a slider input for number of genera to show 
    sidebarLayout(
        sidebarPanel(

            # Input: Select a file ----
            fileInput("file1", "Choose CSV File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
            # Input: Select which barcode to show
            # depends on output.barcode
            h4("Barcode"),
            uiOutput("selectBarcode"),

            # Input: Slider for # of genera to put in plot
            sliderInput("reads_cutoff",
                        "Minimum number of reads:",
                        min = 1,
                        max = 50,
                        value = 30),
            
            # Horizontal line ----
            tags$hr(style="border-color: black;"),
            h4("csv file info:"),
            # Input: Checkbox if file has header ----
            checkboxInput("header", "File has header", TRUE),
        
            # Input: Select separator ----
            radioButtons("sep", "Type of Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
            
            # Input: show input file?
            checkboxInput("showInput", "Show input file", FALSE),
            
             # Horizontal line ----
            tags$hr(style="border-color: black;"),
            h4("Table display:"),
            
            # Input: Select number of rows to display ----
            #h5("Note: only select All if file is small."),
            radioButtons("disp", "Rows to Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head"),
            
            # Ouput: download Excel file of the table or png/pdf of the plot
            h5("Download Data"),
            downloadButton("downloadTable", "Download Table"),
            
            radioButtons("fileOutType", label = "Select plot file type", choices=list("png", "pdf")),
            downloadButton("downloadPlot", "Download Plot")

        ),
    
        # Main panel for displaying outputs ----
        mainPanel(
            #h2("Genus Counts"),
            # Ouput: histogram ----
            plotOutput("countPlot"),
            tags$hr(style="border-color: black;"),
            # Output: Show input data file ----
            # How to show filename here?
            h4("Input file"),
            textOutput("file1_name"),
            tableOutput("input_csv"),
            tags$hr(style="border-color: black;"),
            # Output: Counts table ----
            h4("Output table"),
            tableOutput("counts_table")
        )
    )
)

# Define server logic required to plot the histogram and show tables
server <- function(input, output) {
    # Allow larger files to be uploaded
    options(shiny.maxRequestSize=50*1024^2) 
    
    inputTable <- reactive ({
        req(input$file1)
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        tryCatch(
            {
                reads <- read.csv(input$file1$datapath,
                                  header = input$header,
                                  sep = input$sep,
                                  na.strings=c("-1","NA"))
                reads <- na.omit(reads) # get rid of missing barcodes, genera
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        return(reads)
    })
    
    # get the list of barcodes from input file for the dropdown menu
    output$selectBarcode <- renderUI ({
        reads = inputTable()
        selectInput("barcode", "Select:", sort(unique(reads$barcode)))
    })
    
    countTable <- function() {
        reads = inputTable()
        
        # count the reads by barcode and by genus
        # calculate the % of barcode reads
        # get the taxon name
        
        counts <- reads %>%
            select(barcode, genus_taxid) %>%
            count(barcode, genus_taxid, name="genus_reads") %>%
            group_by(barcode) %>%
            mutate(barcode_reads=sum(genus_reads)) %>%
            mutate(genus = getTaxonomy(genus_taxid, sqlFile = "nameNode.sqlite",
                                     desiredTaxa = c("genus"))) 
            # get the species(?)
            #mutate(species)=getTaxonomy(genus_taxid, sqlFile = "nameNode.sqlite",
            #              desiredTaxa = c("species"))

        # get the 16S rRNA copy number
        copynumdb <- get16Sdb() # from getDBs.R
        # taxid	rank	name	childcount	min	max	mode	median	mean	stddev	sum16slist
        # According to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4021573/pdf/2049-2618-2-11.pdf
        # r = c/g x 100/S, where r = adjusted reads, c = counts, g = copy #, S=community richness
        counts <- left_join(counts, copynumdb, by = c("genus_taxid" = "taxid"))
        counts$mean[is.na(counts$mean)] <- 1 # change the NAs to 1
        # 8/12/20 divide reads by copy number, use that to calculate % of barcode
        counts <- counts %>%
            mutate(genus_reads_norm = genus_reads/mean) %>%
            mutate(perc_of_barcode = 100*genus_reads_norm/barcode_reads)

        # get pathogen data here
        # genusTaxid Genus Pathogen Reference
        # Pathogen = SP - severe pathogen, P - pathogen, B - biofilm, N - nonpathogen
        pathDB <- getPathogenDB() # from getDBs.R
        # reads in as text, convert genus_taxid to numeric for join
        pathDB <- transform(pathDB, genus_taxid=suppressWarnings(as.numeric(genus_taxid)))
        counts <- left_join(counts, pathDB, by="genus_taxid")
        counts$Pathogen[is.na(counts$Pathogen)] <- "U" # make NAs U for Unknown
        return(counts)
    }
    
    createPlot <- function () {
        # wait until we have an input file
        req(input$file1)
        # wait until input$barcode is set
        if(is.null(input$barcode)){return(NULL)}
        
        counts <- countTable()
        # select out the barcode and min # of reads
        plot_counts <- counts %>%
            filter(genus_reads >= input$reads_cutoff & barcode==input$barcode)
        # make the plot object
        group.colors <- c(N = "green", P = "blue", SP ="red", B = "light blue", U = "gray")
        group.labels <- c(N = "Nonpathogen", P = "Pathogen", SP = "Severe pathogen", B = "Biofilm", U = "Unknown")
        
        ggplot(plot_counts, aes(genus, perc_of_barcode, fill=Pathogen)) +
            geom_bar(stat="identity") +
            theme(axis.text.x = element_text(angle = 90),
                  legend.text = element_text(size = 8),
                  legend.key.size = unit(4,"mm")) +
            xlab("Genus") +
            ylab("% of Barcode Reads") +
            ggtitle(paste(input$barcode, " (min reads: ", input$reads_cutoff,
                          ")", sep='')) +
            guides(fill = guide_legend(ncol=1),
                   shape = guide_legend(override.aes = list(size = 0.5))) +
            #Specify colours
            scale_fill_manual(values=group.colors, labels=group.labels)
        
            # Method for adding different plot features based on inputs:
            # df<-data.frame(x=1:10,y=11:20)
            # Switch=F
            # ggplot(df,aes(x,y))+
            #{if(Switch)geom_hline(yintercept=15)}+
            #geom_point()
            # You cannot use + inside the {, but you can put everything into a list 
            # if you have multiple steps. 
            # {if(Switch)list(geom_hline(yintercept=15), geom_hline(yintercept = 13))}
    }
    
    output$countPlot <- renderPlot({
        createPlot()
    })
    
    # print out the name of the input file
    output$file1_name <- renderText({ 
       paste("file name: ", input$file1$name)
    })
   
    # 8/12/20 - displaying the input file is not working...
    output$input_csv <- renderDataTable({

        if (input$showInput == TRUE) {
            reads = inputTable()
            return(head(reads))
        }
    },
    # attempt to format column widths
    options = list(
        autoWidth = TRUE,
        columnDefs = list(list(width = '200px', targets = "_all")))
    )
    # for selected columns: targets = c(1,3)
    # To set different column widths for multiple columns you can use: 
    #columnDefs = (list(list(width = '200px', targets =c(0, 2)), 
    #list(width = '80px', targets =c(6))))
    
    outputTable <- function() {
        if(is.null(input$barcode)){return(NULL)}
        countsOut <- countTable()
        # need to remove some of the columns here
        #barcode	genus_taxid	genus_reads	barcode_reads	genus	rank	name	childcount	min	max	mode	median	mean	stddev	sum16slist	genus_reads_norm	perc_of_barcode	Genus	Pathogen	
        countsOut <- select(countsOut, barcode, genus_taxid, 'Genus name' = genus, Pathogen, 
                            genus_reads, barcode_reads, 'Normalized genus reads' = genus_reads_norm, 
                            '% of barcode' = perc_of_barcode )
        countsOut <- countsOut[countsOut$barcode==input$barcode,]
        countsOut <- countsOut[order(-countsOut$genus_reads),]
        return(countsOut)
    }
    
    output$counts_table <- renderTable({
        # display the output table
        countsOut = outputTable()
        if(input$disp == "head") {
            return(head(countsOut))
        }
        else {
            return(countsOut)
        }
        
    })
    
    ### Download the plot
    output$downloadPlot <- downloadHandler(
        filename = function(){
            paste(input$barcode, input$fileOutType, sep='.')
        },
        content = function(file) {
            if (input$fileOutType=="png")
                png(file)
            else
                pdf(file)
            print(createPlot())
            dev.off()
    })
    
    # Download Excel version of output table ----
    output$downloadTable <- downloadHandler(
        filename = function() {
            paste(input$barcode, ".xlsx", sep = "")
        },
        content = function(file) {
            countsOut <- outputTable()
            write_xlsx(countsOut, file)
        }    
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
