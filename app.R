library(bslib);
library(png);
library(dplyr);
#library(arrow);
library(tidyr);
library(shiny);
library(shinyjs);
library(ggplot2);
library(gridExtra);
library(data.table)
library(DT);
library(kableExtra);
library("formattable");
library(UpSetR);
library(ComplexHeatmap);

# Add cosmo theme and customise for AGRF with CSS
agrf <- bs_theme(version=5, bootswatch="spacelab", tab_font_size = "36px") %>%   #bootswatch="cosmo"
  bs_add_rules(list(sass::sass_file("www/bootstrap.scss")))



# Define UI #####################################################
ui <- fluidPage(
  theme = agrf,
  useShinyjs(),
  
  titlePanel(title = "", windowTitle =  " Find my SNP "),
  
  # Add AGRF banner at the top
  headerPanel(
    div(
      style = "display: flex; align-items: center; justify-content: space-between;",
      div(
        style = "text-align: left; border: 1px solid #ccc; padding: 10px; border-radius: 5px;",
        p(
          style = "font-size: 12px;", # Adjust font size here
          strong("Contact our Customer Care Team:"),
          br(),
          "Freecall 1300 247 301",
          br(),
          "Phone +61 3 9321 3700",
          br(),
          a("CustomerCare@agrf.org.au", href = "mailto:CustomerCare@agrf.org.au"),
          br(),
          a("www.agrf.org.au", href = "https://www.agrf.org.au")
        )
      ),
      a(
        img(src = "agrflogo.png", width = "200px", height = "auto"), # Adjust logo size here
        href = "https://agrf.org.au", target = "_blank"
      )
    )
  ),
  
  
  navbarPage(
    tags$style(HTML("
          table.table-striped {
                width: 100% !important;
                margin-left: 0 !important;   /* Align table to the left */
                padding: 0 !important;       /* Remove any padding */
                margin-right: 0 !important;  /* Remove any right margin */                
            }
         ")),
    
    
    title = tags$span(imageOutput("logo", height = "10", inline = FALSE), " - ", style = "color: #eeeeee; font-size: 2px; font-weight: regular;"),
    
    tabPanel(HTML("<span style='font-size:24px; color:#009FDF;'>Genes</span>"),
             
             sidebarLayout(
               sidebarPanel(
                 style = "height: 550px;",
                 h6("Type or paste refSeq gene id(s) into the field separated by commas. Note: the search is case sensitive."),
                 textInput("geneInput", "Enter Gene(s):", ""),
                 actionButton("plotButton", "Generate plots"),
                 br(),
                 br(),          
                 tableOutput("smalltable"),
                 br(),
                 
                 uiOutput("dynamicPlotHeight"),
                 plotOutput("barsPlots", height = 100)
               ),
               mainPanel(
                 column(width = 12, align = "left", 
                        uiOutput("venDiagramText"),
                        plotOutput("venDiagram", width = "100%", height = "450px") 
                 ),
                 h6("The EPIC array numbers on the UsetR diagram correspond to nearby SNPs, not the actual methylation probes on the EPIC array. Those SNPs are located within a distance of about 1kb around the methylation probes. The number of methylation probes for the EPIC array within the queried gene(s) presented in the bar chart figures on the left side."),
                 br(),
                 hr(),
                 h6("Table with all matched SNPs detected in all selected genes"),
                 DTOutput("longtable")
               ),
               fluid = FALSE
             )),
    tabPanel(HTML("<span style='font-size:24px; color:#009FDF;'>SNPs</span>"),
             sidebarLayout(
               sidebarPanel(
                 textAreaInput("manualInput", "Enter SNPs ID manually", ""),
                 hr(),
                 fileInput("file", "Upload SNP file (one rs-id per row)"),
                 uiOutput("file_summary"),
                 br(),
                 
                 actionButton("submitSNPs", "Submit"),
                 #               h6(style = "font-style: italic; color: #c8c8c8;", "Note: This process may take several minutes to search through over 13 million SNP records. Please be patient."),
                 br(),
                 hr(),
                 
                 plotOutput("heatmap_plot", height = 600)
               ),
               mainPanel(
                 h5(a("UpSetR  diagram", href = "https://upset.app/")),
                 column(width = 12, align = "left", 
                        plotOutput("snpsupset", width = "100%", height = "450px") 
                 ),
                 h6("The EPIC array number on the UsetR diagram and the heatmap correspond to nearby SNPs, not the actual methylation probes on the EPIC array."),
                 br(),
                 hr(),
                 h6("Table with all matched SNPs detected in all selected genes"),
                 DTOutput("snpstable"),
                 
               )
             )
    ),
    tabPanel(HTML("<span style='font-size:24px; color:#009FDF;'>Info</span>"),  
             titlePanel(tags$h4("MethylationEPIC v2 and selected genotyping content search application"), windowTitle = "AGRF app info"),
             mainPanel(
               column(width = 12, align = "left", 
                      HTML("
              <p>This R-Shiny application provides a user-friendly interface for searching, comparing, and analysing the content of Illumina's mainstream genotyping and methylation chips. This application searches for Gene and SNP IDs. It consists of three main tabs: 'Genes', 'SNPs', and 'Info'.</p>
              
              <p><strong>Genes Tab:</strong> Users can input RefSeq gene IDs to generate bar plots showing the number of associated markers: SNPs for genotyping arrays and CpG sites for the EPIC methylation v2.0 chip. The search looks for letter patterns; case-sensitive multiple matches may occur, especially for genes with numerical suffixes. Detailed information about markers with matched gene symbols is displayed in a table. The UpSet plot shows the total number of matched SNPs for the EPIC array. The methylation array does not use probes for SNPs, they are not physically present on the EPIC chip. The SNPs shown are RS-IDs collected within a 1Kb distance for all CpG probes on the array, matched for gene symbols. As such, you may observe more hits in your UpSet plot than shown on the bar chart on the left side for the EPIC array. This approach for loose SNPs search was undertaken in an effort to determine the number of SNPs located near EPIC's methylation probes that intersect with the genotyping chip for the gene of interest. If such shared SNPs exist, the intersect number will be shown on the vertical bars with rows at the bottom of the figure connected by vertical lines.</p>
              
              <p><strong>SNPs Tab:</strong> Users can input SNP IDs manually or upload a file at this tab. The application searches for exact matches for these SNPs across genotyping chip content and SNPs associated with EPIC's CpG probes. Matched annotation is displayed in a table. A heatmap visualisation shows marker counts across genes (rows) and chips (columns). The UpSet plot on this tab shows only exact matches for the SNPs list entered or uploaded by the user. The search is performed across over 13M RS-IDs associated with the EPIC chip and can be slow.</p>
              
              <p><strong>Info Tab:</strong> Presents a summary table with information about two microarray content, including gene coverage and additional details. 

              <p>Overall, the application simplifies comparisons of microarray contents by providing visualisations and detailed statistics for particular selections of genes and SNPs of interest.</p>
              
          "),
                      br(),
                      
                      tableOutput("chipinfolong_table"), 
                      
                      br(),
                      hr(),
                      HTML("
              <p><strong>Note:</strong> The chip annotation used for chromosomes and coordinates is <strong>hg19</strong>. Array annotation files were downloaded from Illumina website and were modified to make uniform dataset. </p>
              <p> The array files were harmonized across seven genotyping chips by using rs-ids as a key. This means that gene names and transcripts from seven illumina genotyping chip's were collected, combined, and re-incorporated back to the chip records for each rs-id. </p>
              <p> The Illumina sourced EPIC array annotation was only modified for the Chr and MapInfo fields. The original chromosome coordinate numbers were converted from hg38 to hg19 using R libraries: rtracklayer and GenomicRanges, with a downloaded chain file from <a href='https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz'  target='_blank'> UCSC </a>.   <p> Shiny application code and additional files and documentation available on <a href='https://github.com/RustTurakulov/GTBUNDLE' target='_blank'> GitHUB </a></p>
              <hr></hr>
              <p> <strong> Disclaimer: </strong></p>
              <p> <i>This Shiny application is in the development stage and is provided for generalized informational and demonstration purposes. The Australian Genome Research Facility (AGRF) and the developer of this application accept no responsibility or liability for any decisions made based on the results or information provided by this application. Users should independently verify the accuracy and reliability of the data presented and exercise their own judgment when interpreting the results. The use of this application is at the user's own risk.</i><p>
          ")
               )
             )
    )  
  )
);

# Create reactive expressions
barPlots           <- reactiveVal(NULL)
all_SEARCHRESULTS  <- reactiveVal(NULL)
SNPSRESULTS        <- reactiveVal(NULL)


# Define server logic
server <- function(input, output) {
  shinyjs::useShinyjs();
  

  output$logo <- renderImage({
    list(src = "agrflogo.png",  width = "220px", height = "50px", alt = "logo")
  }, deleteFile = FALSE);
  
  
  #  output$chipinfolong_table <- renderDT({
  output$chipinfolong_table <- function() {
    req(chipinfolong)
    # Convert Chip column to HTML links
    chipinfolong$Chip <- sprintf("<a href='%s' target='_blank'>%s</a>", 
                                 c("https://sapac.illumina.com/products/by-type/microarray-kits/infinium-methylation-epic.html",
                                   "https://sapac.illumina.com/products/by-type/microarray-kits/infinium-global-screening.html",
                                   "https://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/infinium-global-screening-array-data-sheet-370-2016-016.pdf",
                                   "https://www.illumina.com/products/by-type/microarray-kits/infinium-psycharray.html",
                                   "https://www.illumina.com/products/by-type/microarray-kits/infinium-asian-screening.html",
                                   "https://www.illumina.com/products/by-type/clinical-research-products/infinium-cytosnp-850k.html",
                                   "https://sapac.illumina.com/products/by-type/microarray-kits/infinium-global-diversity.html",
                                   "https://sapac.illumina.com/products/by-type/clinical-research-products/infinium-global-diversity-array-cytogenetics-8.html"
                                 ),
                                 c("EPIC", "GSAv3", "GSAv3 MD", "GSAv3 PS", "GSAv3 AS", "CYTO", "GDAv1", "GDAv1 CYTO"));
    
    chipinfolong$Genes_hg19 <- color_bar("#ED8800")(as.numeric(chipinfolong$Genes_hg19))
    chipinfolong$Sites     <- color_bar("#009FDF")(chipinfolong$Sites)
    kbl(chipinfolong, escape = F) %>%
      kable_paper("hover", full_width = T) %>%
      column_spec(1, width = "2cm") %>%
      column_spec(2, width = "3cm") %>%
      column_spec(3, width = "3cm") %>%
      column_spec(4, width = "3cm") %>%
      column_spec(5, width = "8cm") %>%
      kable_styling("striped", full_width = TRUE, bootstrap_options = c("striped",  "hover", "condensed"), position = "left") %>%    #"bordered",
      as.character()      
  }
  
  
  # Function to generate SNP search results table
  generateSNPTable <- function(searchResults) {
    if (!is.null(searchResults) && length(searchResults) > 0) {
      flattened_RESULTS <- do.call(rbind, lapply(searchResults, as.data.frame))
      snp_table <- data.frame(SNP_ID = generate_unique_rs_ids(flattened_RESULTS$Name), Chip = flattened_RESULTS$Chip)
      return(snp_table)
    } else {
      return(data.frame(SNP_ID = character(), Chip = character()))
    }
  }
  
  
  output$file_summary <- renderUI({
    req(input$file)
    file <- input$file
    if (is.null(file))
      return(NULL)
    
    df <- read.csv(file$datapath, header = TRUE)
    
    # Check if the file contains at least one column
    if (ncol(df) < 1) {
      return(NULL)
    }
    # Assuming the first column contains SNP IDs
    snp_ids <- df[[1]]
    snp_counts <- length(unique(snp_ids));
    summary_table <- data.frame( loaded = "unique snp ids", nn = snp_counts)
    renderTable(summary_table); ## Just a check:  snps list is loaded 
    
  })
  
  
  observeEvent(input$submitSNPs,  {
    # Combine SNP IDs from manual input and file upload
    manual_snps <- if (input$manualInput != "") unlist(strsplit(input$manualInput, "\\s+")) else character(0)
    file_snps <- if (!is.null(input$file)) {
      df <- read.csv(input$file$datapath, header = FALSE)
      if (ncol(df) >= 1) unique(df[[1]]) else character(0)
    } else {
      character(0)
    }
    snp_ids       <- unique(c(manual_snps, file_snps))
    snp_counts    <- length(unique(snp_ids));
    summary_table <- data.frame( detected = "unique snp ids", nn = snp_counts )
    
    
    # Search for SNPs
    snptablesextract <- do.call(rbind, search_snps_table(hg19, snp_ids));  
    SNPSRESULTS=c();
    SNPSRESULTS$EPIC      <- unlist(snptablesextract[which(snptablesextract$Chip == "epic"),]$Name);    
    SNPSRESULTS$GSAv3     <- unlist(snptablesextract[which(snptablesextract$Chip == "gsav3"),]$Name);    
    SNPSRESULTS$GSAv3md   <- unlist(snptablesextract[which(snptablesextract$Chip == "gsav3md"),]$Name);    
    SNPSRESULTS$GSAv3ps   <- unlist(snptablesextract[which(snptablesextract$Chip == "gsav3ps"),]$Name);        
    SNPSRESULTS$GSAv3as   <- unlist(snptablesextract[which(snptablesextract$Chip == "gsav3as"),]$Name);
    SNPSRESULTS$CYTO      <- unlist(snptablesextract[which(snptablesextract$Chip == "cyto"),]$Name);
    SNPSRESULTS$GDAv1     <- unlist(snptablesextract[which(snptablesextract$Chip == "gdav1"),]$Name);    
    SNPSRESULTS$GDAv1cyto <- unlist(snptablesextract[which(snptablesextract$Chip == "gdav1cyto"),]$Name);    
    
    output$snpsupset <- renderPlot({
      req(input$submitSNPs)
      upset(fromList(SNPSRESULTS), main.bar.color = "#ED8800", sets.bar.color="navy", nsets = 8,  text.scale=2)
    })
    
    
    ### snpstable
    output$snpstable <- renderDT({
      if ( length(snp_ids) > 0) {
        req(snptablesextract);
        datatable(snptablesextract,
                  caption = ' Long list of search results across annotation of two chips in vertical format. ',
                  extensions = 'Buttons',
                  options = list(autoWidth = FALSE,
                                 lengthMenu = list(c(10, 50, 100, 500, -1), c('10', '50', '100', '500', 'Full view for export')),
                                 scrollY = TRUE,
                                 scrollX = TRUE,
                                 columnDefs = list(list(visible = FALSE, targets = c(0,1,3))),
                                 paging = TRUE,
                                 searching = TRUE,
                                 fixedColumns = TRUE,
                                 ordering = TRUE,
                                 dom = 'Blfrtip',
                                 pageLength = 10, ## length of the table
                                 buttons = list(
                                   list(extend = 'colvis', text = 'Column Visibility'),
                                   list(extend = 'collection', text = 'Export', 
                                     buttons = list(
                                       list(extend = 'copy',   filename = 'AGRF_findmysnp_hg19table'),
                                       list(extend = 'csv',    filename = 'AGRF_findmysnp_hg19table'),
                                       list(extend = 'excel',  filename = 'AGRF_findmysnp_hg19table'),
                                       list(extend = 'pdf',    filename = 'AGRF_findmysnp_hg19table')
                                  ))
                                 )),
                  rownames  = FALSE) %>%
          formatStyle(columns = c( 1,2 ), fontSize = '75%')
      } else {
        datatable(data.frame(), caption = 'No data available.')
      }
    })
    
    ### heatmap_plot ########################
    snptablesextract <- do.call(rbind, search_snps_table(hg19, snp_ids));
    snptablesextract$`Gene(s)` <- as.character(snptablesextract$`Gene(s)`)
    snptablesextract$Chip      <- as.character(snptablesextract$Chip)
    HEATMX <-  as.data.frame.matrix(table(snptablesextract$`Gene(s)`, snptablesextract$Chip))
    
    if (ncol(HEATMX) < 8) {
      missing_cols <- setdiff(c( "epic", "gsav3", "gsav3md", "gsav3ps", "gsav3as", "cyto", "gdav1", "gdav1cyto"), colnames(HEATMX))
      for (mycol in missing_cols) {
        HEATMX[, mycol] <- rep(0, nrow(HEATMX))
      }
    }
    
    habar = HeatmapAnnotation(
      "SNPs" = anno_barplot(
        round(colSums(HEATMX),1), 
        bar_width = 1, 
        gp = gpar(col = "white", fill = "yellow"), 
        border = TRUE,
        add_numbers = TRUE,
        height = unit(3, "cm")), 
      show_annotation_name = TRUE,
      annotation_name_rot  = 0,
      annotation_name_side = "left"
    )
    habar@anno_list$"SNPs"@fun@var_env$gp$fill <- "#009FDF"
    custom_palette <- colorRampPalette(c("white", "gray", "#009FDF"));
    col_fun <- custom_palette(9);
    #    col_fun <- RColorBrewer::brewer.pal(9,"PiYG") #"YlGnBu" "Blues"
    
    ht = Heatmap(as.matrix(HEATMX), 
                 name = "SNPs",
                 col =  col_fun,
                 cluster_columns   = FALSE,
                 show_row_dend     = FALSE, 
                 show_column_dend  = FALSE,
                 show_column_names = TRUE,
                 show_row_names    = TRUE,
                 heatmap_width     = unit(12, "cm"),
                 heatmap_height    = unit(18, "cm"),
                 top_annotation    = c( habar ),
                 column_names_rot  = 45, 
                 column_names_side = "top",
                 heatmap_legend_param = list(direction = "vertical"),
                 column_title_gp = gpar(fontsize = 24),
                 row_names_gp    = gpar(fontsize = 12),
                 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text( HEATMX[i, j], x, y, gp = gpar(fontsize = 14, col="navy"))
                 }
    );
    
    
    output$heatmap_plot <- renderPlot({
      draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left")
    });
    ########################################## end of heatmap
  })     
  
  observeEvent(input$plotButton, {
    genes <- unique(unlist(strsplit(input$geneInput, "\\s*,\\s*")))
    
    # Clear individual plots and tables before generating new 
    individual_plots  <- list();
    all_SEARCHRESULTS(NULL);
    all_search_results <- if (!is.null(all_SEARCHRESULTS())) all_SEARCHRESULTS() else list()
    
    for (gene in genes) {
      results <- generateBarCharts(gene)
      if (!is.null(results[[1]]) && inherits(results[[1]], "gg")) {
        individual_plots  <- c(individual_plots,  list(results[[1]]))
        all_search_results <- c(all_search_results, list(results[[2]]))
      }
    }
    
    barPlots(individual_plots)
    all_SEARCHRESULTS(all_search_results)
  })
  
  output$barsPlots <- renderPlot({
    if (!is.null(barPlots())) {
      # Calculate the height dynamically based on the number of genes
      plot_height <- 200 + length(barPlots()) * 150  # Adjust the multiplier as needed
      
      # Create a list of heights for each plot
      plot_heights <- rep(1, length(barPlots())) * (plot_height / sum(plot_height))
      
      # Print each individual plot in a single column with dynamically adjusted height
      print(grid.arrange(grobs = barPlots(), ncol = 1, heights = plot_heights))
    } else {
      # If no plots, create an empty plot with helix banner
      #      plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
      bg_image <- readPNG("www/agrfwebbanner.png")
      plot_height <- 150;
      scatter_plot <- ggplot(data = data.frame(x=c(0,1), y=c(1,0)), aes(x = x, y = y)) +
        geom_point(shape = 21, colour = "white", fill = "white", size = 0.2) +  
        annotation_custom(rasterGrob(bg_image, width = unit(1, "npc"), height = unit(1, "npc")), 
                          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
        theme_void() +  # Remove default theme elements
        theme(plot.margin = margin(0, 0, 0, 0))
      print(scatter_plot)
    }
  }, height = function() {
    # Set the overall height dynamically based on the number of genes
    200 + length(barPlots()) * 150  # Adjust the multiplier as needed
  });
  
  
  
  output$longtable <- renderDT({
    if (!is.null(all_SEARCHRESULTS()) && length(all_SEARCHRESULTS()) > 0) {
      flattened_SEARCHRESULTS <- do.call(rbind, lapply(all_SEARCHRESULTS(), as.data.frame))
      datatable(flattened_SEARCHRESULTS,
                caption = ' long list of search results across chips annotations in vertical format ',
                extensions = 'Buttons',
                options = list(autoWidth = FALSE,
                               lengthMenu = list(c(10, 50, 100, 500, -1), c('10', '50', '100', '500', 'Full view for export')),
                               scrollY = TRUE,
                               scrollX = TRUE,
                               columnDefs = list(list(visible = FALSE, targets = c(1:3))),
                               paging = TRUE,
                               searching = TRUE,
                               fixedColumns = TRUE,
                               ordering = TRUE,
                               dom = 'Blfrtip',
                               pageLength = 10,
                               buttons = list(
                                 list(extend = 'colvis', text = 'Column Visibility'),
                                 list(extend = 'collection', text = 'Export', 
                                      buttons = list(
                                        list(extend = 'copy',   filename = 'AGRF_findmysnp_hg19table'),
                                        list(extend = 'csv',    filename = 'AGRF_findmysnp_hg19table'),
                                        list(extend = 'excel',  filename = 'AGRF_findmysnp_hg19table'),
                                        list(extend = 'pdf',    filename = 'AGRF_findmysnp_hg19table')
                                      ))
                               )),
                rownames  = FALSE) %>%
        formatStyle(columns = c(2, 4), fontSize = '75%')
    } else {
      datatable(data.frame(), caption = 'No data available.')
    }
  })
  
  output$venDiagramText <- renderUI({
    if (is.null(all_SEARCHRESULTS())) {
      fluidRow(
        column(
          width = 12,
          h4("Use this search application to understand the best genotyping array for your candidate genes."),
          h5("Test Instructions:"),
          tags$ol(
            tags$li("Enter the symbol or multiple symbols of the gene(s) you wish to search for in the ", tags$code("Enter Gene(s):"), " text input box on the left side of the application. For example" , tags$code("NF1, BRAF") ),
            tags$li("Click the ", tags$strong("Generate plots"), " button located below the input box.")
          ),
          h5("Expected Output:"),
          tags$ol(
            tags$li(
              tags$strong("Bar Plots:"),
              tags$ul(
                tags$li("Bar plots will display the count of markers associated with the entered gene(s) on each microarray chip (GDA or EPIC)."),
                tags$li("Note: on the barplots GDA array shows unique rs-ids matched to gene symbol and  for EPIC array bar shows number of cpg probes for the gene.")
              )
            ),
            tags$li(
              tags$strong("UpSetR Diagram:"),
              tags$ul(
                tags$li("The UpSetR diagram illustrates the overlap of markers between the EPIC and GDA arrays."),
                tags$li("This diagram is rs id oriented no cpg probes from EPIC array used. RS-id collected from EPIC annotation within 1kb distances from all gene matched cpg probes.")
              )
            ),
            tags$li(
              tags$strong("Table with Matched SNPs:"),
              tags$ul(
                tags$li("A table will show detailed information about the matched markers, including their ID, chip, chromosome, and map info."),
                tags$li("You can export table or copy it to the buffer, however the only visible rows and columns will be collected."),
                tags$li("Table column and row number visibility can be adjusted at the left side corner of the table.") 
              )
            )
          ),
          h5("Contacts and credits:"),
          tags$ol(
            tags$strong("AGRF Sales:"),
            tags$ul(
              tags$li("Desley Pitcher: desley.pitcher@agrf.org.au")
            ),             
            tags$strong("AGRF Genotyping:"),
            tags$ul(
              tags$li("Melinda Ziino: melinda.ziino@agrf.org.au")
            ),
            tags$strong("Developers:"),
            tags$ul(
              tags$li("Lesley Gray: lesley.gray@agrf.org.au"),
              tags$li("Rust Turakulov: rust.turakulov@agrf.org.au")
            )
          )
        )
      )
    }else{
      h5(a("UpSetR  diagram", href = "https://upset.app/"))
      
    }
  })
  
  output$venDiagram <- renderPlot({
    if (!is.null(all_SEARCHRESULTS())) {
      venn_table <- do.call(rbind, lapply(all_SEARCHRESULTS(), as.data.frame));
      venn_epic  <- venn_table[which(venn_table$Chip == "EPIC"),]
      venn_data  <- list(
        EPIC     = generate_unique_rs_ids(venn_epic$Name),
        GSA      = as.vector(venn_table[which(venn_table$Chip == "GSA"),     "Name"]),
        GSA_md   = as.vector(venn_table[which(venn_table$Chip == "GSA_md"),  "Name"]),
        GSA_ps   = as.vector(venn_table[which(venn_table$Chip == "GSA_ps"),  "Name"]),
        GSA_as   = as.vector(venn_table[which(venn_table$Chip == "GSA_as"),  "Name"]),
        CYTO     = as.vector(venn_table[which(venn_table$Chip == "CYTO"),    "Name"]),
        GDA      = as.vector(venn_table[which(venn_table$Chip == "GDA"),     "Name"]),
        GDA_cyto = as.vector(venn_table[which(venn_table$Chip == "GDA_cyto"),"Name"])
      )
      
      
      upset(fromList(venn_data), main.bar.color =  "#ED8800", sets.bar.color="navy", nsets = 8,  text.scale=2)  # "#009FDF"
      
    } else {
      plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
    }
  })
  
  output$smalltable <- renderTable({
    flattened_SEARCHRESULTS <- if (!is.null(all_SEARCHRESULTS())) {
      do.call(rbind, lapply(all_SEARCHRESULTS(), as.data.frame))
    } else {
      data.frame()
    }
    
    if (nrow(flattened_SEARCHRESULTS) > 0) {
      summary_table_data <- data.frame(Chip = flattened_SEARCHRESULTS$Chip, SNPs = 1)
      summary_table_data <- aggregate(SNPs ~ Chip, data = summary_table_data, length)
      colnames(summary_table_data) <- c("Chip", "Probes")  
      
      summary_table_data <- merge(summary_table_data, chipinfolong[, c("Chip", "Name")], by="Chip", all.x=TRUE)
      
      return(summary_table_data[order(-summary_table_data$Probes), ])
    } else {
      return(data.frame(Chip = character(), Probes = integer()))
    }
  });
  
}






#Trouble with slow AWS load
options(
  shiny.server = list(
    app_init_timeout = 120,  # 60 seconds is default
    app_idle_timeout = 300   # 3600 seconds (1 hour)
  )
)

hg19           <- c();
#epic_csv <- open_dataset(sources = "epic.csv", format = "csv");
#epic_csv       <- open_dataset("epic_parquet/"); #faster
#hg19$epic      <- as.data.frame(epic_csv);

hg19$epic      <- fread("epic.csv");
hg19$gsav3     <- fread("gsav3.csv");
hg19$gsav3md   <- fread("gsav3md.csv");
hg19$gsav3ps   <- fread("gsav3ps.csv");
hg19$gsav3as   <- fread("gsav3as.csv");
hg19$cyto      <- fread("cyto.csv");
hg19$gdav1     <- fread("gdav1.csv");
hg19$gdav1cyto <- fread("gdav1cyto.csv");

EPIC      <- setDT(hg19$epic);
GSAv3     <- setDT(hg19$gsav3);
GSAv3md   <- setDT(hg19$gsav3md);
GSAv3ps   <- setDT(hg19$gsav3md);
GSAv3as   <- setDT(hg19$gsav3as);
CYTO      <- setDT(hg19$cyto);
GDAv1     <- setDT(hg19$gdav1);
GDAv1cyto <- setDT(hg19$gdav1cyto);


generate_unique_rs_ids <- function(column) {
  concatenated_rs_ids  <- paste(column, collapse = ",")
  individual_rs_ids    <- unlist(strsplit(concatenated_rs_ids, ","))
  unique_rs_ids        <- as.character(factor(individual_rs_ids, levels = unique(individual_rs_ids)))
  return(unique_rs_ids)
};

chipinfolong = data.frame(
  "Chip" = c(
    "EPIC",
    "GSA",
    "GSA_md",    
    "GSA_ps",
    "GSA_as",
    "CYTO",
    "GDA",
    "GDA_cyto"),
  
  "Name" = c(
    " MethylationEPIC v2.0 ",
    " Global Screening Array v3 ",
    " Global Screening Array v3 + MultiDisease",
    " Global Screening Array v3 + Psychology",
    " Global Screening Array v3 + Asian panel",  
    " Cyto850K V1.4",
    " Global Diversity Array v1",
    " Global Diversity Array with Cytogenetics" ),
  
  "Genes_hg19" = c( 26447,  29310,   30076,   30076,  37954,  39725,   52593,   54715 ),
  "Sites"    =  c( 935000, 654027,  654027,  654027, 659184, 971836, 1904599, 1977761 ),
  
  "Info"= c(
    "Robust methylation profiling microarray with extensive coverage of CpG islands, genes, and enhancers. Use for epigenome-wide association studies. There are 13,086,857 known SNPs in vicinity of methylation sites probed on the array.",
    "Includes a multiethnic genome-wide backbone, clinical research variants, QC markers, and custom add-on content.",
    "This next-generation genotyping array provides a scalable, cost-effective solution for population genetics, pharmacogenomics studies, and precision medicine research.",    
    "Additional ~50K genetic markers associated with common psychiatric disorders. Focus on multiple diseases and selected by the Psychiatric Genomics Consortium from The Broad Institute in Boston, USA, have been incorporated into the core GSA array.",
    "The Infinium Asian Screening Array (ASA) combines genome-wide coverage of East Asian populations.A powerful, cost-effective genotyping array for large-scale genetic studies and pharmacogenomics in East Asian populations. Fixed markers: ~ 660,000, Custom marker add-on capacity up to 50,000",
    "This consortium-built array provides comprehensive coverage of cytogenetically relevant genes for congenital disorders and cancer research",
    "The Global Diversity Array (GDA) BeadChip combines exceptional coverage of clinical research variants with optimized multi-ethnic, genome-wide content. Approximately 1.8M markers",
    "Profile clinical research variants associated with congenital disease, cancer research, pharmacogenomics (PGx), and exome content with this microarray. Approximately 1.8M markers"
  )
);      


# Function to generate bar plots
generateBarCharts <- function(gene) {
  EPIClong  <- EPIC[grep(paste0("^", gene, "|,\\s+", gene),   EPIC$"Gene(s)", perl = TRUE), ]
  EPICshort <- EPIClong %>%
    group_by(`Transcript(s)`) %>%
    summarise(
      Chr             = first(Chr),
      MapInfo         = first(MapInfo),      
      Name            = paste(Name, collapse = ","),
      "Transcript(s)" = first(`Transcript(s)`),    
      "Gene(s)"       = first(`Gene(s)`)
#      misc            = first(misc)
    ) %>%
    ungroup()
  
  
  SEARCHRESULTS <- list(
    EPIC     = EPICshort,
    GSA      = GSAv3[  grep(paste0("^", gene, "|,\\s+", gene),    GSAv3$"Gene(s)",     perl = TRUE), ],
    GSA_md   = GSAv3md[grep(paste0("^", gene, "|,\\s+", gene),    GSAv3md$"Gene(s)",   perl = TRUE), ],    
    GSA_ps   = GSAv3ps[grep(paste0("^", gene, "|,\\s+", gene),    GSAv3ps$"Gene(s)",   perl = TRUE), ],
    GSA_as   = GSAv3as[grep(paste0("^", gene, "|,\\s+", gene),    GSAv3as$"Gene(s)",   perl = TRUE), ],
    CYTO     = CYTO[   grep(paste0("^", gene, "|,\\s+", gene),    CYTO$"Gene(s)",      perl = TRUE), ],
    GDA      = GDAv1[  grep(paste0("^", gene, "|,\\s+", gene),    GDAv1$"Gene(s)",     perl = TRUE), ],
    GDA_cyto = GDAv1cyto[grep(paste0("^", gene, "|,\\s+", gene),  GDAv1cyto$"Gene(s)", perl = TRUE), ]
  )
  
  chip_data <- data.frame(cbind(Chip = names(SEARCHRESULTS), SNP_count = as.numeric(lapply(SEARCHRESULTS, nrow))))
  chip_order <- c( "EPIC", "GSA", "GSA_md", "GSA_ps", "GSA_as", "CYTO", "GDA", "GDA_cyto");
  chip_data$Chip <- factor(chip_data$Chip, levels = chip_order)
  chip_data$SNP_count <- as.numeric(chip_data$SNP_count)
  # Find the maximum value for adjusting the y-axis range
  max_value <- max(chip_data$SNP_count)
  
  bar_plots <- ggplot(chip_data, aes(x = Chip, y = SNP_count, fill = Chip)) +
    geom_bar(stat = "identity", fill = "#009FDF") +
    geom_label(aes(label = SNP_count), vjust = -0.2, fontface = "bold", color = "navy", size = 5,
               fill = "#ED8800", label.padding = unit(0.4, "lines")) +
    labs(title = paste(gene, " - gene has markers / probes"), x = " ", y = "SNP Count") +
#    scale_fill_manual(values = c( "EPIC" = "orange",          "GSA" = "deepskyblue1", "GSA_ps"="palegoldenrod", "GSA_md" = "azure2",
#                                "GSA_as" = "lavenderblush",  "CYTO" = "cornsilk",     "GDA"="paleturquoise3",  "GDA_cyto"= "wheat" )) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +  # Adjust the expand values as needed
    coord_cartesian(ylim = c(0, max_value * 1.2)) +  # Adjust the multiplier as needed
    theme_minimal() +
    theme(axis.title.y = element_text(face = "bold")) +
    theme(axis.text.x = element_text(face = "bold", size = 14)) +
    guides(fill = "none")
  
  # Collecting some info for the summary data
  non_empty_tables <- sapply(SEARCHRESULTS, function(x) nrow(x) > 0)
  
  if (any(non_empty_tables)) {
    # At least one table is non-empty
    selected_data_frames <- lapply(names(SEARCHRESULTS)[non_empty_tables], function(chip) {
      xtable <- SEARCHRESULTS[[chip]]
      xtable$Chip <- chip
      return(xtable)
    })
    
    # Bind the rows of the data frames into a single data frame
    tabledata <- do.call(rbind, selected_data_frames)
  } else {
    # All tables are empty, create an empty data frame with the same column structure
    tabledata <- data.frame(Chip = character(), stringsAsFactors = FALSE)
  }
  
  tabledata$MapInfo <- as.numeric(tabledata$MapInfo)
  return(list(bar_plots = bar_plots, chip_data = tabledata))
};


search_snps_table <- function(hg19, snp_ids) {
  search_results <- list()
  if (length(snp_ids) == 0) {
    return(search_results)
  }
  for (table_name in names(hg19)) {
    table <- hg19[[table_name]]
    if ("Name" %in% colnames(table)) {
      matches <- table[table$Name %in% snp_ids, ]
      if (nrow(matches) > 0) {  # Check if matches are found
        matches$Chip <- table_name
        search_results[[table_name]] <- matches
      }
    }
  }
  return(search_results)
}




# Run the application
shinyApp(ui, server)