#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Michael Sierk
# 7/29/20
# Last update: 1/6/20
#   changed app to SequaPath, reconfigured diretories, redeployed to shinyapps.io
#   -> problem during redeployment when clicking "publish" button, used this command at the console:
#       deployApp(appDir=getwd(), appFileManifest="appFileManifest.txt")
#
# For Nephros, Inc.
# 
# Read in .csv from EPI2ME software, count reads from each taxid. Output:
#   1) Table of counts normalized by copy number
#   2) Histogram of percentage of reads in the barcode, labeled by pathogenicity

# csv headers
# filename	read_id	exit_status	runid	taxid	barcode	accuracy	lca	genus_taxid
# Allow user to select file from filesystem
library(shiny)
library(dplyr)
library(ggplot2)
library(taxonomizr) # converts taxid to genus & species
library(writexl) # for writing Excel output
source("getDBs.R") # script for updating DBs
library(shinyBS) # for mouseover popups
library(cowplot)

# Change to TRUE to (re)generate the the taxonomy database nameNode.sqlite
make_taxonomy_db = FALSE
if (make_taxonomy_db) {
    makeTaxonomyDB()
}


# Define UI for application 
#ui <- fluidPage(theme = shinytheme("flatly"),
ui <- fluidPage(
                                
    # Application title
    title = "SequaPath",
    
    fluidRow(
        column(4, height=100, background = "light-blue",
               # Use imageOutput to place the image on the page
               imageOutput("logo", height="50px")
        ),
        column(8, height=100, background = "light-blue",
               h2("SequaPath Bacterial Census")
        )
    ),
    
    fluidRow(
        # Main panel for displaying outputs ----
        column(12,
            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(type = "tabs",
                    tabPanel("Plot", plotOutput("countPlot")),
                    tabPanel("Output Table", htmlOutput("table_barcode"), 
                             tableOutput("counts_table")),
                    tabPanel("Input File", h4("Input file"),
                             textOutput("file1_name"),
                             tableOutput("input_csv"))
                    
            ),
            tags$hr(style="border-color: black;")
        ), # main panel
    ),
    
    fluidRow(
        column(4,
            # Input: Select a file ----
            fileInput("file1", "Choose CSV File (200 MB limit)",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
           # Input: Select which barcode to show
            # depends on output.barcode
            h4("Barcode:"),
            uiOutput("selectBarcode"),
            # Add: need free text entry for experiment # (ID?), and notes
            # - will be associated with Barcode, included in output table
            textInput("barcodeLabel", "Barcode Label", value = "", width = NULL, placeholder = NULL),
           
            # Input: Slider for # of genera to put in plot
            sliderInput("reads_cutoff",
                        "Minimum number of reads:",
                        min = 1,
                        max = 50,
                        value = 30),
            bsTooltip("file1", "Input .csv file must have following headers: ",
                     "top", options = list(container = "body"))
        ),
            
        column(4,
               # Horizontal line ----
               h4("Output Table Display"),
            # Input: Select number of rows to display ----
            radioButtons("disp", "Output Rows:",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head", inline=TRUE),
            hr(),
            h4("Input .csv file info"),
            # Input: Checkbox if file has header ----
            checkboxInput("header", "File has header", TRUE),
            
            # Input: Select separator ----
            radioButtons("sep", "Type of Separator:",
                         choices = c(Comma = ",",
                                     Semicolon = ";",
                                     Tab = "\t"),
                         selected = ",", inline = TRUE),
            
            # Input: show input file?
            checkboxInput("showInput", "Show input file head", FALSE)
        ),
        column(4, 
            # Horizontal line ----
            h4("Download options"),
            # Download the report
            div(style = "padding: 4px 4px", 
                downloadButton("report", "Generate report")
            ),
            # Ouput: download Excel file of the table
            #h5("Download Data"),
            div(style = "padding: 4px 4px", 
                downloadButton("downloadTable", "Download Table (Excel)")
            ),   
            # Download the histogram plot
            div(style = "padding: 4px 4px", 
                downloadButton("downloadPlot", "Download Plot"),
                radioButtons("fileOutType", label = "Select plot file type:", 
                            choices=list("png", "pdf"), inline = TRUE)
            )
        ),
    ) # fluidRow
) # ui

# Define server logic required to plot the histogram and show tables
server <- function(input, output, session) {
    # Allow larger files to be uploaded
    options(shiny.maxRequestSize=200*1024^2) 
    
    inputTable <- reactive ({
        req(input$file1)
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to throw error
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
    
    load16Sdb <- reactive ({
        # get the 16S rRNA copy number
        
        copynumdb <- get16Sdb() # from getDBs.R
        return(copynumdb)
    })
    
    loadPathogenDB <- reactive ({
        pathDB <- getPathogenDB(redo_table = FALSE) # from getDBs.R
        return(pathDB)
    })
    
    # add a tooltip 
    addTooltip(session, "file1", "Input .csv file must have following headers: ", 
               placement = "top", trigger = "hover",
               options = NULL)
    
    # get the list of barcodes from input file for the dropdown menu
    output$selectBarcode <- renderUI ({
        reads = inputTable()
        selectInput("barcode", NULL, sort(unique(reads$barcode)))
    })
    
    countTable <- reactive ({
        reads = inputTable()
        # input table header: filename,read_id,exit_status,runid,taxid,barcode,accuracy,lca,genus_taxid
        #   keep taxid, barcode, accuracy, genus_taxid
        #   reads below threshold have taxid, genus_taxid = -1
        # - count the reads by barcode and by genus
        # - divide by copy number if have it
        # - calculate the % of barcode reads
        # - get the taxon name
        
        counts <- reads %>%
            select(barcode, taxid, genus_taxid, accuracy) %>%
            group_by(barcode) %>%
            count(barcode, genus_taxid, name="genus_reads") %>%
            mutate(barcode_reads=sum(genus_reads)) %>%
            mutate(genus = getTaxonomy(genus_taxid, sqlFile = "db/nameNode.sqlite",
                                     desiredTaxa = c("genus"))) # getTaxonomy from taxonomizr
            # get the species(?)
            #mutate(species)=getTaxonomy(genus_taxid, sqlFile = "nameNode.sqlite",
            #              desiredTaxa = c("species"))

        # get the 16S rRNA copy number
        copynumdb <- load16Sdb()
        # taxid	rank	name	childcount	min	max	mode	median	mean	stddev	sum16slist
        # According to PMC4021573
        # r = c/g x 100/S, where r = adjusted reads, c = counts, g = copy #, S=community richness
        counts <- left_join(counts, copynumdb, by = c("genus_taxid" = "taxid"))
        counts$mean[is.na(counts$mean)] <- 1 # change the ones not in the DB to 1
        # 8/12/20 divide reads by copy number, use that to calculate % of barcode
        counts <- counts %>%
            mutate(genus_reads_norm = genus_reads/mean) %>%
            mutate(perc_of_barcode = 100*genus_reads_norm/barcode_reads)

        # get pathogen data here
        # genusTaxid Genus Species Pathogen(N/P1/P2/P3) Biofilm(B) Gram (Pos/Neg) Notes and Synonyms
        # Pathogen = P1 - severe pathogen, P2 - pathogen, P3 = ?, N - nonpathogen
        pathDB <- loadPathogenDB() # from getDBs.R
        #print("pathDB: ")
        #print(head(pathDB))
        # reads in as text, convert genus_taxid to numeric for join
        # Note 8/25/20: genus_taxid is being converted to a float
        pathDB <- transform(pathDB, genus_taxid=suppressWarnings(as.integer(genusTaxid)))
        counts <- left_join(counts, pathDB, by="genus_taxid")
        counts$Pathogen[is.na(counts$Pathogen)] <- "U" # make NAs U for Unknown
        #print("Counts table: ")
        #print(head(counts))
        return(counts)
    })
    
    # select out just the data from selected barcode for plotting/output
    barcodeTable <- reactive ({
        if(is.null(input$barcode)){return(NULL)}
        counts <- countTable()
        #head(counts)
        # full table headers:
        #barcode genus_taxid genus_reads barcode_reads genus rank name childcount min max mode median mean stddev sum16slist genus_reads_norm perc_of_barcode Genus Pathogen	
        countsOut <- counts %>%
            filter(genus_reads >= input$reads_cutoff & barcode==input$barcode) %>%
            select(barcode, genus, Pathogen, Biofilm, 
                    genus_reads, perc_of_barcode, 
                    mean, barcode_reads, 
                    genus_reads_norm) %>%
                    arrange(desc(perc_of_barcode))
         
        return(countsOut)
    })
    
    barcodeReads <- reactive ({
        # want to keep barcode reads as a separate variable
        countsOut <- barcodeTable()
        get_barcode_reads <- countsOut %>% slice(1)
        barcodeReads <- get_barcode_reads$barcode_reads
        return(barcodeReads)
    })
    
    createPlot <- reactive ({
        # wait until we have an input file
        req(input$file1)
        # wait until input$barcode is set
        if(is.null(input$barcode)){return(NULL)}
        
        plot_counts <- barcodeTable()
        
        # make the plot object
        group.colors <- c(P1 ="red", P2 = "blue", P3 = "light blue", N = "light green", U = "gray")
        group.labels <- c(P1 = "CDC Pathogen", P2 = "Nephros Pathogen", P3 = "Koch Institute Pathogen", 
                          N = "Nonpathogen", U = "Unknown")
        # col = Biofilm, plus scale_color_manual outlines genera that are Biofilm producers
        pl <- ggplot(plot_counts, aes(genus, perc_of_barcode, fill=Pathogen, group=Biofilm)) +
            geom_bar(aes(col = Biofilm), stat="identity") +
            theme_cowplot(12) +
            theme(plot.title = element_text(size = 18, face = "bold"),
                  plot.subtitle = element_text(size = 14),
                  plot.caption = element_text(size = 14, face = "italic"),
                  axis.text.x = element_text(angle = 90, size = 14, hjust = 1, vjust = 0.5),
                  axis.title.x = element_text(face = "bold", size = 18),
                  axis.text.y = element_text(size = 14),
                  axis.title.y = element_text(face = "bold", size = 18),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 14, face = "bold"),
                  legend.key.size = unit(6,"mm")) +
            xlab("Genus") +
            ylab(paste0("% of Barcode Reads\n(min reads: ", 
                        input$reads_cutoff,") ")) +
            guides(fill = guide_legend(ncol=1),
                   shape = guide_legend(override.aes = list(size = 0.5))) +
            #Specify colours
            scale_fill_manual(values=group.colors, labels=group.labels) +
            scale_color_manual(values=c(B="black"), guide = "none") +
            labs(title = input$barcode, subtitle = input$barcodeLabel, caption = "Outline = Biofilm producer")
            #annotate("text", Inf, Inf, label = "Outline = Biofilm producer", hjust = 1, vjust = 1, size = 5)
        
        pl <- ggdraw(pl) + 
                draw_image("nephros_logo.jpg", x = 0.95, y = 0.1, hjust = 1, vjust = 1, width = 0.1, height = 0.1)
        return(pl)
    })
    
    # display the logo with renderImage
    output$logo <- renderImage({
        # Return a list containing the filename and alt text
        list(src = "nephros_logo.jpg", alt = "logo", height=50)
    }, deleteFile = FALSE)
    
    output$countPlot <- renderPlot({
        pl <- createPlot()
        print(pl)
    })
    
    # print out the name of the input file
    output$file1_name <- renderText({ 
       paste("file name: ", input$file1$name)
    })
   
    # display the input table
    output$input_csv <- renderTable({

        if (input$showInput == TRUE) {
            reads = inputTable()
            return(head(reads))
        }
    },
    # attempt to format column widths
    #options = list(
    #    autoWidth = TRUE,
    #    columnDefs = list(list(width = '200px', targets = "_all")))
    )
    # for selected columns: targets = c(1,3)
    # To set different column widths for multiple columns you can use: 
    #columnDefs = (list(list(width = '200px', targets =c(0, 2)), 
    #list(width = '80px', targets =c(6))))
    
    # print out the barcode and label before rendering the table
    output$table_barcode <- renderUI({ 
        barcodeReads <- barcodeReads()
        str1 <- paste(input$barcode,"(",input$barcodeLabel,")")
        str2 <- paste("Barcode reads: ", barcodeReads)
        HTML(paste(str1, str2, sep='<br/>'))
    })
    
    output$counts_table <- renderTable({
        # display the output table
        table_counts = barcodeTable()
        get_barcode_reads <- distinct(table_counts, barcode_reads)
        barcodeReads <- get_barcode_reads$barcode_reads
        print(paste("Barcode reads: ", barcodeReads))
        table_counts <- table_counts %>%
                        ungroup() %>%
                        select('Genus' = genus, Pathogen, Biofilm, '% of barcode' = perc_of_barcode,
                        'Normalized genus reads' = genus_reads_norm, genus_reads, 
                        'Copy number' = mean)
        if(input$disp == "head") {
            return(head(table_counts))
        }
        else {
            return(table_counts)
        }
        
    })
    
    ### Print a report
    output$report <- downloadHandler(
        filename = "report.pdf",
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path(tempdir(), "report.Rmd")
            file.copy("report.Rmd", tempReport, overwrite = TRUE)
            tempLogo <- file.path(tempdir(), "nephros_logo.png")
            file.copy("nephros_logo.png", tempLogo, overwrite = TRUE)
            tempLogo <- file.path(tempdir(), "nephros_logo.jpg")
            file.copy("nephros_logo.jpg", tempLogo, overwrite = TRUE)
            tempsty <- file.path(tempdir(), "unicode-math.sty")
            file.copy("unicode-math.sty", tempsty, overwrite = TRUE)
            
            # Set up parameters to pass to Rmd document
            # Need to figure out what parameters to send here
            # the plot 
            # the input file name
            # the barcode (multiple plots for different barcodes?)
            # table output for top 10 genera?
            table_counts = barcodeTable()
            barcodeReads = barcodeReads()
            table_counts <- table_counts %>%
                ungroup() %>%
                select('Genus' = genus, Pathogen, Biofilm, '% of barcode' = perc_of_barcode,
                       'Normalized genus reads' = genus_reads_norm, 'Genus reads' = genus_reads, 
                       'Copy number' = mean) %>% 
                        mutate_if(is.numeric, format, digits=2)
            
            params <- list(pl = createPlot(), inpName = input$file1$name, barcode=input$barcode,
                           barcodeLabel = input$barcodeLabel, barcodeReads = barcodeReads, 
                           table_counts = table_counts)
            
            # Knit the document, passing in the `params` list, and eval it in a
            # child of the global environment (this isolates the code in the document
            # from the code in this app).
            rmarkdown::render(tempReport, output_file = file,
                              params = params,
                              envir = new.env(parent = globalenv())
            )
        }
    )

    ### Download the plot
    output$downloadPlot <- downloadHandler(
        filename = function(){
            paste(input$barcode, input$fileOutType, sep='.')
        },
        content = function(file) {
            if (input$fileOutType=="png")
                png(file, width=10, height=7, units="in", res=300)
            else
                pdf(file, width=10, height=7)
            pl <- createPlot()
            print(pl)
            dev.off()
    })
    
    # Download Excel version of output table ----
    output$downloadTable <- downloadHandler(
        filename = function() {
            paste(input$barcode, ".xlsx", sep = "")
        },
        content = function(file) {
            table_counts <- barcodeTable()
            table_counts <- table_counts %>%
                ungroup() %>%
                select('Genus' = genus, Pathogen, Biofilm, '% of barcode' = perc_of_barcode,
                       'Normalized genus reads' = genus_reads_norm, 'Genus reads' = genus_reads, 
                       'Copy number' = mean, 'Barcode reads' = barcode_reads)
            write_xlsx(table_counts, file)
        }    
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
