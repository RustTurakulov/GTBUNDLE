library(dplyr);
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


hg19gda  <- readRDS("hg19_GDA.rds");
#ahg19epic <- readRDS("hg19_EPIC.rds");
hg19       <- hg19gda; 
hg19$epic  <- readRDS("hg19_EPIC.long.rds")
#hg19 <- c(hg19gda, hg19epic)

GDA  <- setDT(hg19$gdav1);
EPIC <- setDT(hg19$epic);


generate_unique_rs_ids <- function(column) {
  concatenated_rs_ids  <- paste(column, collapse = ",")
  individual_rs_ids    <- unlist(strsplit(concatenated_rs_ids, ","))
  unique_rs_ids        <- as.character(factor(individual_rs_ids, levels = unique(individual_rs_ids)))
  return(unique_rs_ids)
};
#epicuniverse <- generate_unique_rs_ids(hg19$epic$Name); 


chipinfolong = data.frame(
  "Chip" = c(
        "GDA",
        "EPIC"),

    "Name" = c(
        "Infinium Global Diversity Array v1",
        "Infinium MethylationEPIC v2.0"),

  "Genes_hg19" = c(32394,  26447 ),
  
  "Sites" = c(1904599, 935000),

  "Info"= c(
       "The Global Diversity Array (GDA) BeadChip combines exceptional coverage of clinical research variants with optimized multi-ethnic, genome-wide content. Approximately 1.8M markers",
       "Robust methylation profiling microarray with extensive coverage of CpG islands, genes, and enhancers. Use for epigenome-wide association studies. There are 13,086,857 known SNPs in vicinity of methylation sites probed on the array."
    )
);      


# Function to generate bar plots
generateBarCharts <- function(gene) {
  EPIClong <- EPIC[grep(paste0("^", gene, "|,\\s+", gene),   EPIC$"Gene(s)", perl = TRUE), ]
  EPICshort <- EPIClong %>%
    group_by(`Transcript(s)`) %>%
    summarise(
      Chr             = first(Chr),
      MapInfo         = first(MapInfo),      
      Name            = paste(Name, collapse = ","),
      "Transcript(s)" = first(`Transcript(s)`),    
      "Gene(s)"       = first(`Gene(s)`),
      misc            = first(misc)
    ) %>%
    ungroup()
  
  
  SEARCHRESULTS <- list(
    GDA      = GDA[grep(paste0("^", gene, "|,\\s+", gene),    GDA$"Gene(s)", perl = TRUE), ],
    EPIC     = EPICshort
  )
  
  chip_data <- data.frame(cbind(Chip = names(SEARCHRESULTS), SNP_count = as.numeric(lapply(SEARCHRESULTS, nrow))))
  chip_order <- c("GDA", "EPIC" )
  chip_data$Chip <- factor(chip_data$Chip, levels = chip_order)
  chip_data$SNP_count <- as.numeric(chip_data$SNP_count)
  # Find the maximum value for adjusting the y-axis range
  max_value <- max(chip_data$SNP_count)
  
  bar_plots <- ggplot(chip_data, aes(x = Chip, y = SNP_count, fill = Chip)) +
    geom_bar(stat = "identity") +
    geom_label(aes(label = SNP_count), vjust = -0.2, fontface = "bold", color = "yellow", size = 6,
               fill = "gray", label.padding = unit(0.2, "lines")) +
    labs(title = paste(gene, " - gene has markers / probes"), x = " ", y = "SNP Count") +
    scale_fill_manual(values = c( "GDA" = "orange", "EPIC" = "deepskyblue1" )) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +  # Adjust the expand values as needed
    coord_cartesian(ylim = c(0, max_value * 1.2)) +  # Adjust the multiplier as needed
    theme_minimal() +
    theme(axis.title.y = element_text(face = "bold")) +
    theme(axis.text.x = element_text(face = "bold", size = 12)) +
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
  
  # Filter GDA data
  gda_matches <- hg19$gdav1[hg19$gdav1$Name %in% snp_ids, ]
  if (nrow(gda_matches) > 0) {  
    gda_matches$Chip <- "GDA"
    search_results$GDA <- gda_matches
  }
  
  # Filter EPIC data
  epic_matches <-  hg19$epic[hg19$epic$Name %in% snp_ids, ]
  if (nrow(epic_matches) > 0) {  
    epic_matches$Chip <- "EPIC"
    search_results$EPIC <- epic_matches
  }
  
  return(search_results)
}





# Define UI #####################################################
ui <- navbarPage(
    tags$style(HTML("
           table.table-striped {
                width: 100% !important;
                margin-left: 0 !important; /* Align table to the left */
                padding: 0 !important; /* Remove any padding */
                margin-right: 0 !important; /* Remove any right margin */                
            }
         ")),

  
  title = tags$span(imageOutput("logo", height = "10", inline = FALSE), " GT CHIPS ", style = "color: white; font-size: 2px; font-weight: regular;"),
  
  tabPanel("Genes",
#           imageOutput("logo", height = "10", inline = TRUE),
           sidebarLayout(
             sidebarPanel(
               style = "height: 275px;",
               h5("Enter refSeq gene id(s) into the field separated by commas. Note the search is case sensitive."),
               textInput("geneInput", "Enter Gene(s):", ""),
               actionButton("plotButton", "Generate plots"),
               
               br(),
               tableOutput("smalltable"),
               hr(),
               
               uiOutput("dynamicPlotHeight"),
               plotOutput("barsPlots", height = 300)
             ),
             mainPanel(
               h4(a("UpSetR  diagram", href = "https://upset.app/")),
               column(width = 12, align = "left", 
                      uiOutput("venDiagramText"),
                      plotOutput("venDiagram", width = "100%", height = "450px") 
               ),
               
               br(),
               hr(),
               h4("Table with all matched SNPs detected in all selected genes"),
               DTOutput("longtable")
             ),
             fluid = FALSE
           )),
  tabPanel("SNPs",
           sidebarLayout(
             sidebarPanel(
               textAreaInput("manualInput", "Enter SNPs ID manually", ""),
               hr(),
               fileInput("file", "Upload SNP file (one rs-id per row)"),
               uiOutput("file_summary"),
               br(),

               actionButton("submitSNPs", "Submit"),
               h5(style = "font-style: italic; color: #555;", "Note: This process may take several minutes to search through over 13 million SNP records. Please be patient."),
               
               br(),
               hr(),
               
               
               plotOutput("heatmap_plot", height = 600)
             ),
             mainPanel(
               h4(a("UpSetR  diagram", href = "https://upset.app/")),
               column(width = 12, align = "left", 
                      plotOutput("snpsupset", width = "100%", height = "450px") 
               ),
               
               br(),
               hr(),
               h4("Table with all matched SNPs detected in all selected genes"),
               DTOutput("snpstable"),

             )
           )
  ),
  tabPanel("Info", 
       titlePanel(tags$h4("MethylationEPIC v2 and Global Diversity Array v1 content search application"), windowTitle = "GT CHIPS"),
       mainPanel(
         column(width = 12, align = "left", 
         HTML("
              <p>The Shiny application provides a user-friendly interface for searching, comparing, and analyzing the content of Illumina's mainstream genotyping and methylation chips. Aplication focused on searching for Gene and SNP IDs. It consists of three main tabs: 'Genes', 'SNPs', and 'Info'.</p>
              
              <p><strong>Genes Tab:</strong> Users can input RefSeq gene IDs to generate bar plots showing the number of associated markers: SNPs for GDA array and CPg sites for EPIC chip. The search looks for letter patterns, and case-sensitive multiple matches may occur, especially for genes with numerical suffixes. Detailed information about markers with matched gene symbol is displayed in a table. Please not be confused with the upset plot which shows total number of matched SNPs for EPIC array. As methylation array do not normally probes SNPs those are not physically present on EPIC chip. Those are just rs-ids collected within 1Kb distance for all cpg probes on array matched for gene symbol. That is the reasoon why you may observe more hits for upsed plot than shown on barchart for the EPIC array. This done for the attemp to find number of SNPs nearby EPIC's methylation probes intersects with GDA chip. If those shared SNPs exist they number will be shown on vertical bar with  connected vertical line.</p>
              
              <p><strong>SNPs Tab:</strong> Users can input SNP IDs manually or upload a file at this tab. The application searches for exact matches for these SNPs across GDA content and SNPs associated with EPIC's cpg probes. Matched annotation displayed in a table. A heatmap visualization shows marker counts across genes (rows) and chips (columns). Ther upset plot on this tab show only exact match for the SNPs list enered by user or uploaded. The serch is performed acros over 13M or rs-ids associated with EPIC chips and can be slow.</p>
              
              <p><strong>Info Tab:</strong> Presents a summary table with information about two microarrays content, including gene coverage and additional details.</p>
              
              <p>Overall, the application simplifies comparisons of  microarray contents by providing visualizations and detailed statistics for particular selection of gene and snps of interest.</p>
              
          "),
         br(),
         
         tableOutput("chipinfolong_table"), 
         
         br(),
         hr(),
         HTML("
              <p><strong>Note:</strong> The chip annotation used for chromosomes and coordinates is <strong>hg19</strong>. Both GDA and EPIC array annotation files were downloaded from Illumina website. </p>
              <p> The GDA array file was harmonized across seven genotyping chips by using rs-ids as keys. This means that gene names and transcripts from seven illumina genotyping chip's were collected, combined, and re-incorporated back to GDA chip records for each rs-id. </p>
              <p> The Illumina sourced EPIC array annotation was only modified for the Chr and MapInfo fields. The original chromosome coordinate numbers were converted from hg38 to hg19 using R libraries: rtracklayer and GenomicRanges, with a downloaded chain file from <a href='https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz'  target='_blank'> UCSC </a>.   <p> Shiny application code and additional files and documentation available on <a href='https://github.com/RustTurakulov/GTBUNDLE' target='_blank'> GitHUB </a></p>
              <hr></hr>
              <p> <strong> Disclaimer: </strong></p>
              <p> <i>This Shiny application is in the development stage and is provided for generalized informational and demonstration purposes. The Australian Genome Research Facility (AGRF) and the developer of this application accept no responsibility or liability for any decisions made based on the results or information provided by this application. Users should independently verify the accuracy and reliability of the data presented and exercise their own judgment when interpreting the results. The use of this application is at the user's own risk.</i><p>
          ")
         )
       )
))

# Create reactive expressions
barPlots           <- reactiveVal(NULL)
all_SEARCHRESULTS  <- reactiveVal(NULL)
SNPSRESULTS        <- reactiveVal(NULL)


# Define server logic
server <- function(input, output) {
  shinyjs::useShinyjs()

  output$logo <- renderImage({
    list(src = "agrflogo.png",  width = "220px", height = "50px", alt = "logo")
  }, deleteFile = FALSE);

  output$showSpinner <- renderUI({
    if (showSpinner()) {
      div(class = "loader")
    } else {
      NULL
    }
  })
  

#  output$chipinfolong_table <- renderDT({
   output$chipinfolong_table <- function() {
      req(chipinfolong)
     # Convert Chip column to HTML links
     chipinfolong$Chip <- sprintf("<a href='%s' target='_blank'>%s</a>", 
                                  c("https://sapac.illumina.com/products/by-type/microarray-kits/infinium-global-diversity.html",
                                    "https://sapac.illumina.com/products/by-type/microarray-kits/infinium-methylation-epic.html"
                                    ),
                                  c("GDA", "EPIC"));
     
      chipinfolong$Genes_hg19 <- color_bar("orange")(as.numeric(chipinfolong$Genes_hg19))
      chipinfolong$Sites     <- color_bar("deepskyblue")(chipinfolong$Sites)
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
    snptablesextract <- do.call(rbind, search_snps_table(hg19, snp_ids));  ## Slow function need to fix 
    SNPSRESULTS=c();
    SNPSRESULTS$GDA  <- unlist(snptablesextract[which(snptablesextract$Chip == "GDA"),]$Name);
    ALLEPICRSID      <- generate_unique_rs_ids(snptablesextract[which(snptablesextract$Chip == "EPIC"), ]$Name);
    SNPSRESULTS$EPIC <- ifelse(length(snp_ids[snp_ids %in% ALLEPICRSID]) > 0, NA, snp_ids[snp_ids %in% ALLEPICRSID])

    output$snpsupset <- renderPlot({
      req(input$submitSNPs)
      upset(fromList(SNPSRESULTS), main.bar.color = "royalblue", sets.bar.color="navy", nsets = 2,  text.scale=2)
    })
    
    
    ### snpstable
    output$snpstable <- renderDT({
      if ( length(snp_ids) > 0) {
        req(snptablesextract);
        datatable(snptablesextract,
                  caption = ' Long list of search results across annotation of two chip in vertical format ',
                  extensions = 'Buttons',
                  options = list(autoWidth = FALSE,
                                 lengthMenu = list(c(10, 50, 100, 500, -1), c('10', '50', '100', '500', 'All')),
                                 scrollY = TRUE,
                                 scrollX = TRUE,
                                 columnDefs = list(list(visible = FALSE, targets = c(0,3,4))),
                                 paging = TRUE,
                                 searching = TRUE,
                                 fixedColumns = TRUE,
                                 ordering = TRUE,
                                 dom = 'Blfrtip',
                                 pageLength = -1, 
                                 buttons = list(
                                    list(extend = 'colvis', text = 'Column Visibility'),
                                    list(extend = 'copy', text = 'Copy',   filename = 'AGRF_GTChips_matchtable_hg19'),
                                    list(extend = 'csv', text = 'CSV',     filename = 'AGRF_GTChips_matchtable_hg19'),
                                    list(extend = 'excel', text = 'Excel', filename = 'AGRF_GTChips_matchtable_hg19')
                                  #  list(extend = 'pdf', text = 'PDF',     filename = 'AGRF_GTChips_matchtable_hg19')
                                )),
                                rownames  = FALSE) %>%
                         formatStyle(columns = c(2, 4), fontSize = '75%')
      } else {
        datatable(data.frame(), caption = 'No data available.')
      }
    })
    
    ### heatmap_plot ########################
    snptablesextract <- do.call(rbind, search_snps_table(hg19, snp_ids));
    snptablesextract$`Gene(s)` <- as.character(snptablesextract$`Gene(s)`)
    snptablesextract$Chip      <- as.character(snptablesextract$Chip)
    HEATMX <-  as.data.frame.matrix(table(snptablesextract$`Gene(s)`, snptablesextract$Chip))
    if (ncol(HEATMX) < 2) {
      missing_cols <- setdiff(c("GDA", "EPIC"), colnames(HEATMX))
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
    habar@anno_list$"SNPs"@fun@var_env$gp$fill <- "royalblue"
    col_fun <- RColorBrewer::brewer.pal(9,"PiYG") #"YlGnBu" "Blues"

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
                 column_names_rot  = 0, 
                 column_names_side = "top",
                 heatmap_legend_param = list(direction = "vertical"),
                 column_title_gp = gpar(fontsize = 20),
                 row_names_gp    = gpar(fontsize = 10),

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
      plot_height <- 300 + length(barPlots()) * 150  # Adjust the multiplier as needed
      
      # Create a list of heights for each plot
      plot_heights <- rep(1, length(barPlots())) * (plot_height / sum(plot_height))
      
      # Print each individual plot in a single column with dynamically adjusted height
      print(grid.arrange(grobs = barPlots(), ncol = 1, heights = plot_heights))
    } else {
      # If no plots, create an empty plot
      plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, ann = FALSE)
    }
  }, height = function() {
    # Set the overall height dynamically based on the number of genes
    300 + length(barPlots()) * 150  # Adjust the multiplier as needed
  });
  
  
  output$longtable <- renderDT({
    if (!is.null(all_SEARCHRESULTS()) && length(all_SEARCHRESULTS()) > 0) {
      flattened_SEARCHRESULTS <- do.call(rbind, lapply(all_SEARCHRESULTS(), as.data.frame))
      datatable(flattened_SEARCHRESULTS,
                caption = ' long list of search results across both chips annotations in vertical format ',
                extensions = 'Buttons',
                options = list(autoWidth = FALSE,
                               lengthMenu = list(c(10, 50, 100, 500, -1), c('10', '50', '100', '500', 'All')),
                               scrollY = TRUE,
                               scrollX = TRUE,
                               columnDefs = list(list(visible = FALSE, targets = c(1:3))),
                               paging = TRUE,
                               searching = TRUE,
                               fixedColumns = TRUE,
                               ordering = TRUE,
                               dom = 'Blfrtip',
                               pageLength = -1, 
                               buttons = list(
                                 list(extend = 'colvis', text = 'Column Visibility'),
                                 list(extend = 'copy', text = 'Copy',   filename = 'AGRF_GTBundle_matchtable_hg19'),
                                 list(extend = 'csv', text = 'CSV',     filename = 'AGRF_GTBundle_matchtable_hg19'),
                                 list(extend = 'excel', text = 'Excel', filename = 'AGRF_GTBundle_matchtable_hg19'),
                                 list(extend = 'pdf', text = 'PDF',     filename = 'AGRF_GTBundle_matchtable_hg19')
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
          h4("Test Instructions:"),
          tags$ol(
            tags$li("Enter the symbol or multiple symbols of the gene(s) you wish to search for in the ", tags$code("Enter Gene(s):"), " text input box on the left side of the application. For example" , tags$code("NF1, BRAF") ),
            tags$li("Click the ", tags$strong("Generate plots"), " button located below the input box.")
          ),
          h4("Expected Output:"),
          tags$ol(
            tags$li(
              tags$strong("Bar Plots:"),
              tags$ul(
                tags$li("Bar plots will display the count of markers associated with the entered gene(s) on each microarray chip (GDA or EPIC)."),
                tags$li("Note: on the barplots GDA array shows unique rs-ids matched to gene symbol and  for EPIC array bar show number of cpg probes for the gene).")
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
                tags$li("You can export table or copy it to the buffer however the only visible rows and columns will be collected"),
                tags$li("Table column and row number visibility can be adjusted at the left side corner of the table.") 
              )
            )
         ),
         h4("Contacts and credits:"),
           tags$ol(
                tags$strong("Developers:"),
                tags$ul(
                  tags$li("Rust Turakulov: rust@agrf.org.au "),
                  tags$li("Lesley Gray: lesley.gray@agrf.org.au ")
                ),
                tags$strong("Genotyping team:"),
                  tags$ul(
                   tags$li("Melinda Ziino: melinda.ziino@agrf.org.au ")
                  )
               )
             )
          )
    }
  })
  
  output$venDiagram <- renderPlot({
    if (!is.null(all_SEARCHRESULTS())) {
      venn_table <- do.call(rbind, lapply(all_SEARCHRESULTS(), as.data.frame));
      venn_epic  <- venn_table[which(venn_table$Chip == "EPIC"),]
      venn_data  <- list(
        GDA  = as.vector(venn_table[which(venn_table$Chip == "GDA"),    "Name"]),
        EPIC = generate_unique_rs_ids(venn_epic$Name)
      )
      
      upset(fromList(venn_data), main.bar.color = "royalblue", sets.bar.color="navy", nsets = 2,  text.scale=2)
     

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
      colnames(summary_table_data) <- c("Chip", "SNPs")  
      summary_table_data <- merge(summary_table_data, chipinfolong[, c("Chip", "Name")], by="Chip")
      
      
      return(summary_table_data)
    } else {
      return(data.frame(Chip = character(), SNPs = integer()))
    }
  });

}

# Run the application
shinyApp(ui, server)