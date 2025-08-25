library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(shinyFeedback)
library(DT)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(plotly)
library(shinyjqui)

library(biomaRt)
library(vegan)
library(tidyr)
library(dplyr)
library(uwot) # umap
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(AnnotationDbi) # kegg
library(DESeq2)


# ===== Functions =====
check_id_type <- function(id_vector) {
  ids <- head(id_vector, 5)
  if (all(grepl("^ENSG[0-9]+(\\.[0-9]+)?$", ids))) {
    return("gene_id")
  } else if (all(grepl("^ENST[0-9]+(\\.[0-9]+)?$", ids))) {
    return("transcript_id")
  } else {
    return("unknown")
  }
}


# ===== UI Main =====
ui <- dashboardPage(
  dashboardHeader(title = "KGE-seq",
                  tags$li(
                    class = "dropdown",
                    style = "padding: 1px; margin-right: 4px;",
                    tags$span("Beta 1.0", style = "font-size: 10px; color: black;"))
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("home")),
      menuItem("Functions", tabName = "FUNCTIONS", icon = icon("dna"),
               menuSubItem("Differental analysis", tabName = "DEanalysis", icon = icon("angle-right"))

      )
    )
  ),
  dashboardBody(
    # theme # Black and Grey
    tags$head(tags$style(HTML('
              /* logo */
              .skin-blue .main-header .logo {
                background-color: #3c3c3c;
                color: #ffffff;
                font-weight: bold;
              }
              .skin-blue .main-header .logo:hover {
                background-color: #3c3c3c;
              }
              
              /* navbar */
              .skin-blue .main-header .navbar {
                background-color: #3c3c3c;
              }
              
              /* sidebar */
              .skin-blue .main-sidebar {
                background-color: #4a4a4a;
              }
              
              /* active sidebar item */
              .skin-blue .main-sidebar .sidebar .sidebar-menu .active a {
                background-color: #6e6e6e;
                color: #ffffff;
                font-weight: bold;
              }
              
              /* sidebar item */
              .skin-blue .main-sidebar .sidebar .sidebar-menu a {
                background-color: #4a4a4a;
                color: #ffffff;
              }
              
              /* sidebar item hover */
              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover {
                background-color: #5a5a5a;
                color: #ffffff;
              }
              
              /* sidebar toggle hover */
              .skin-blue .main-header .navbar .sidebar-toggle:hover {
                background-color: #5a5a5a;
              }
              
              /* body content */
              .content-wrapper, .right-side {
                background-color: #fafafa;
                color: #000000;
              }
              
              /* box color */
              .box {
                background-color: #ffffff;
                border: 1px solid #cccccc;
              }
              
              /* box header */
              .box-header {
                background-color: #4a4a4a;
                color: #ffffff;
                border-bottom: 1px solid #cccccc;
                font-weight: bold;
              }
          '))),
    
    
    tabItems(
      # home
      tabItem("Home",
              fluidRow(
                box(title="Search Genes", width=12, solidHeader=TRUE, collapsible=TRUE,
                    textInput("ensembl_id", "Ensembl ID or Symbol:", value="ENSG00000139618", width="20%"),
                    selectInput("org_db_home", "Select organism",
                                choices = c("Human" = "hsapiens_gene_ensembl",
                                            "Mouse" = "mmusculus_gene_ensembl",
                                            "Rat"   = "rnorvegicus_gene_ensembl"),
                                selected = "hsapiens_gene_ensembl", width="20%"),
                    loadingButton("search_gene", "Search gene", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    HTML('<br><br>'),
                    fluidRow(
                      column(width=6, verbatimTextOutput("gene_info"))
                    )
                  ) # box
              ) # Row
      ), # Item
      
      # DE analysis
      tabItem("DEanalysis",
              fluidRow(
                box(title="The reqired format", width=12, solidHeader=TRUE, collapsible=TRUE,
                    helpText("Please ensure the data format match the demo table"),
                    HTML("<br>"),
                    fluidRow(
                      column(6,
                             tabPanel("Metadata",
                                      HTML("<ol>\
                                          <li>Accession_ID : The unique identifier for each sample</li>\
                                          <li>Labels or Patients : A group info assigned to each sample</li>\
                                          </ol>"),
                                      div(DT::dataTableOutput("demo_meta")))),
                      column(6,
                             tabPanel("Expression matrix",
                                      HTML("<ol>\
                                          <li>Row: Accession_ID</li>\
                                          <li>Column: Ensembl gene IDs</li>\
                                          </ol>"),
                                      div(DT::dataTableOutput("demo_matrix"))))
                    )
                  ), # box
                
                # Basic analysis
                box(title="Basic analysis", width=12, solidHeader=TRUE, collapsible=TRUE,
                    # upload meta
                    fluidRow(
                      column(6, 
                             div(style="height:85px", 
                                 fileInput(inputId = "meta_file",label = HTML("Input your metadata with csv format :"), 
                                           accept = c(".csv"), placeholder = "metadata.csv", multiple = F)))),
                    # upload matrix
                    fluidRow(
                      column(6,
                             div(style="height:85px",
                                 fileInput(inputId = "matrix_file",label = HTML("Input your dataset with csv format :"), 
                                           accept = c(".csv"), placeholder = "dataset.csv", multiple = F)))),
                    HTML("<br>"),
                    loadingButton("upload_files", "Upload files", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    HTML("<br><br>"),
                    fluidRow(
                      column(6, uiOutput("given_meta_ui")),
                      column(6, uiOutput("given_matrix_ui"))
                    ),
                    # read depth
                    HTML("<br><br>"),
                    fluidRow(
                      column(8, uiOutput("reads_depth"))
                    ),
                    # PCA
                    HTML("<br><br><br>"),
                    conditionalPanel(
                      condition = "output.reads_depth !== null && output.reads_depth !== undefined",
                      HTML("<h4><b>Distribution of expression</b></h4>"),
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          selectInput("distr_type", "Method",
                                      choices = c("PCA" = "pca", 
                                                  "UMAP" = "umap")),
                          selectInput("distr_group", "Group by",
                                      choices = c()),
                          HTML("<br>"),
                          loadingButton("distr_button", "Get results", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                        ),
                        mainPanel(
                          plotOutput("distr_plot")
                        ))
                      ),
                    HTML("<br><br>"),
                    conditionalPanel(
                      condition = "output.r_depth !== null && output.r_depth !== undefined",
                      downloadButton("export_basic", label = "Export results"),
                      style = "position: absolute; bottom: 10px; right: 10px;"
                    )
                    
                ), # box
                
                # DE analysis
                box(title="Differential analysis", width=12, solidHeader=TRUE, collapsible=TRUE,
                    radioButtons("DE_method", "Differential Methods", c("DESeq2"='deseq'), selected="deseq"),
                    selectInput("org_db_DE", "Select organism", choices=c("Human", "Mouse", "Rat"), width="20%"),
                    conditionalPanel(
                      condition = "input.matrix_file.length > 0",
                      fluidRow(column(6, checkboxGroupInput("DE_group", "Design by", choices = NULL)))
                    ),
                    loadingButton("run_DE", "Run", loadingLabel = "Processing...", style="color: #444; background-color: #f4f4f4; border-color: #ddd;"),
                    HTML("<br><br>"),
                    uiOutput("DE_finish"),
                    HTML("<br><br>"),
                    tabsetPanel(
                      # volcano plot module
                      tabPanel("Enriched genes",
                               sidebarLayout(
                                 sidebarPanel(
                                   width = 3,
                                   numericInput("log2FC", label='log2 fold change cutoff', value=1, min=0, max=50, step=1),
                                   numericInput("q_val", label='q-value cutoff', value=0.05, min=0, max=5, step=1),
                                   selectInput("DE_name", label='Focus on', choices=""),
                                   HTML("<br>"),
                                   loadingButton("vol_button", "Get results", loadingLabel = "Processing...",
                                                 style = "color: #444; background-color: #f4f4f4; border-color: #ddd;")
                                 ),
                                 mainPanel(
                                   uiOutput("top_genes_table_ui")
                                 )
                               ),
                               fluidRow(
                                 column(
                                   width = 12,
                                   plotlyOutput("vol_plot", height = "600px")
                                 )
                               ),
                               HTML("<br><br>"),
                               useShinyjs(),
                               downloadButton("export_volplot", "Export results",
                                              style = "position: absolute; bottom: 10px; right: 10px; z-index: 10;")
                               
                      ),
                      # pathway analysis module
                      tabPanel("Pathway analysis",
                               sidebarLayout(
                                 sidebarPanel(
                                   width = 3,
                                   selectInput("path_db", label="Database", choices=c("GO", "KEGG"), selected="GO"),
                                   uiOutput("dynamic_path_ui"),
                                   HTML("<br>"),
                                   loadingButton("path_button", "Get results", loadingLabel = "Processing...",
                                                  style = "color: #444; background-color: #f4f4f4; border-color: #ddd;")
                                   ),
                                 mainPanel(
                                   uiOutput("top_path_table_ui")
                                 )
                               ),
                               fluidRow(
                                 column(
                                   width = 12,
                                   jqui_resizable(plotOutput("path_plot", height = "600px"))
                                 )
                               ),
                               HTML("<br><br>"),
                               useShinyjs(),
                               downloadButton("export_pathway", "Export results",
                                              style = "position: absolute; bottom: 10px; right: 10px; z-index: 10;")
                      )
                    ) # table set
                ) # box
              ) # Row
      ) # Item
      
      
    ) # items
  ) # body
) # boardpage


# ===== Server Main =====
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  
  # reactive values
  values_upload <- reactiveValues(meta=NULL, mat=NULL)
  values_basic <- reactiveValues(depth_plot=NULL, pca_plot=NULL)
  values_DE <- reactiveValues(method=NULL, fit=NULL, org=NULL, log2FC=NULL, q_val=NULL, res=NULL, vol_plot=NULL, DEG_table=NULL)
  values_path <- reactiveValues(gene_list_ensembl=NULL, gene_list_entrez=NULL, path_dot=NULL, path_table=NULL)
  
  # demo format
  output$demo_meta <- DT::renderDT({
    demoMeta <- data.frame('Accession_ID'=paste0('Sample_',1:4), 
                           'Labels'=c('control', 'CRC', 'control', 'CRC'),
                           'Patients'=c('Patient_A', 'Patient_A', 'Patient_B', 'Patient_B'))
    datatable(demoMeta, options=list(scrollX=T))
  })
  
  output$demo_matrix <- DT::renderDT({
    set.seed(1234)
    genes <- paste0("ENSG00000", sample(100000:999999, 4))
    samples <- paste0("Sample_", 1:4)
    demo_mat <- matrix(sample(0:10, 16, replace = TRUE), nrow = 4, ncol = 4,
                        dimnames = list(genes, samples))
    datatable(as.data.frame(demo_mat), options=list(scrollX=T))
  })
  
  
  # Event 1 search_gene
  observeEvent(input$search_gene, {
    selected_db <- input$org_db_home
    MART <- useMart("ensembl", dataset = selected_db)
    es_id <- trimws(input$ensembl_id)
    
    if (grepl("^ENS", es_id, ignore.case = TRUE)) {
      filter_to_use <- "ensembl_gene_id"
    } else {
      symbol_filter_map <- c(
        "hsapiens_gene_ensembl"    = "hgnc_symbol",
        "mmusculus_gene_ensembl"   = "mgi_symbol",
        "rnorvegicus_gene_ensembl" = "rgd_symbol"
      )
      filter_to_use <- symbol_filter_map[[selected_db]]
    }
    
    gene_info_ <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "chromosome_name", "start_position", "end_position"),
            filters = filter_to_use,
            values = es_id,
            mart = MART)
   
    output$gene_info <- renderText({
      req(gene_info_)
      paste0(gene_info_$external_gene_name, "\n",
             "Ensembl ID: ", gene_info_$ensembl_gene_id, "\n\n",
             "Description: ", sub(" \\[.*", "", gene_info_$description), "\n",
             "Chromosome name: ", gene_info_$chromosome_name, "\n",
             "Start position: ", gene_info_$start_position, "\n",
             "End position: ", gene_info_$end_position, "\n")
    })
    resetLoadingButton("search_gene")
  }) # end of search gene
  
  
  # Event 2 upload_files
  observeEvent(input$upload_files, {
    # check if the files prepare
    if(is.null(input$meta_file) | is.null(input$matrix_file)) {
      showNotification("Uploaded files are not ready.", type = "error")
      resetLoadingButton("upload_files")
      req(FALSE)
    }
    
    # check if it is csv file
    if (!grepl("\\.csv$", input$meta_file$name, ignore.case = TRUE)) {
      showNotification("Metadata file must be a .csv file", type = "error")
      req(FALSE)
      resetLoadingButton("upload_files")
    }
    if (!grepl("\\.csv$", input$matrix_file$name, ignore.case = TRUE)) {
      showNotification("Matrix file must be a .csv file", type = "error")
      req(FALSE)
      resetLoadingButton("upload_files")
    }
    meta <- read.csv(input$meta_file$datapath,
                     check.names = FALSE)
    values_upload$meta <- meta
    raw_mat <- read.csv(input$matrix_file$datapath, 
                        row.names = 1,
                        check.names = FALSE)
    raw_mat <- raw_mat[rowSums(is.na(raw_mat)) == 0, ] # remove invalid gene data
    values_upload$mat <- raw_mat
    
    # check the format
    if (!"Accession_ID" %in% colnames(meta)) {
      showNotification("Metadata must contain 'Accession_ID' column", type = "error")
      req(FALSE)
      resetLoadingButton("upload_files")
    }
    
    if (!any(grepl("^ENSG", rownames(raw_mat)))) {
      showNotification("Row names of matrix must be gene IDs (e.g., starting with 'ENSG')", type = "error")
      req(FALSE)
      resetLoadingButton("upload_files")
    }
    
    # update group selection
    group_info <- setdiff(colnames(meta), "Accession_ID")
    # remove non-str attributes
    group_info[sapply(meta[, group_info, drop = FALSE], function(x) is.character(x) || is.factor(x))]
    updateSelectInput(session, "distr_group", choices=unique(group_info))
    updateCheckboxGroupInput(session, "DE_group", choices=unique(group_info))
    
    # check gene id type
    id_type <- check_id_type(rownames(raw_mat))
    if (id_type == "transcript_id") {
      showNotification("Detect the transcript ID. Please convert it into gene-level.", type="error")
    } else if (id_type == "gene_id") {
      showNotification("Dectect the gene-level matrix.")
    } else {
      showNotification("Cannot detect the type of data.", type="error")
    }
    
    # deal with version number
    gene_ids <- rownames(raw_mat)
    remove_ver <- sub("\\..*", "", gene_ids)
    rownames(raw_mat) <- remove_ver
    
    # 1_basic_output
    mat_plus <- raw_mat + 1
    long_mat <- mat_plus %>% pivot_longer(cols = everything(), names_to = "Sample", values_to = "Counts")
    depth_plot <- ggplot(long_mat, aes(x = Sample, y = Counts)) +
      geom_boxplot(fill="grey") +
      scale_y_continuous(trans = "log2", labels = scales::math_format(2^.x, format = log2)) +
      ylab("log2 (Counts+1)") + theme(axis.text.x = element_text(size = 6.5, angle = 90))
    
    values_basic$depth_plot <- depth_plot
    # Output render
    output$reads_depth <- renderUI({
      if(is.null(depth_plot)) {
        return(NULL)
      }
      
      tagList(
        HTML("<h4><b>Reads depth</b></h4>"),
        plotOutput("r_depth")
      )
    })
    
    output$r_depth <- renderPlot({
      depth_plot
    })
    
    output$given_meta_ui <- renderUI({
      tagList(
        HTML("<h4>Metadata</h4>"),
        div(DT::dataTableOutput("given_meta"))
      )
    })
    
    output$given_matrix_ui <- renderUI({
      tagList(
        HTML("<h4>Expression data</h4>"),
        div(DT::dataTableOutput("given_matrix"))
      )
    })
    
    output$given_meta <- DT::renderDT({
      datatable(meta, options=list(scrollX=T), rownames=FALSE)
    })
    
    output$given_matrix <- DT::renderDT({
      datatable(raw_mat, options=list(scrollX=T))
    })
    
    resetLoadingButton("upload_files")
  }) # end of upload_files

  # 2_PCA_plot
  observeEvent(input$distr_button, {
    selected_method <- input$distr_type
    given_group <- input$distr_group
    meta <- values_upload$meta
    raw_mat <- values_upload$mat
    
    if(selected_method == "pca"){
      log_mat <- log2(raw_mat + 1) # log2(count+1)
      log_mat <- log_mat[apply(log_mat, 1, function(x) var(x) > 0), ]
      pca <- prcomp(t(log_mat), scale. = TRUE)
      pca_df <- data.frame(PC1 = pca$x[,1],
                           PC2 = pca$x[,2],
                           Sample = colnames(log_mat))
      
      pca_df$group <- meta[[given_group]][match(pca_df$Sample, meta$Accession_ID)]
      d_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
                        geom_point(size = 3) + # different shape
                        geom_text_repel(aes(label = Sample), size = 2, show.legend = FALSE) + # label do not overlap
                        theme_minimal() +
                        labs(x = "PC1", y = "PC2") +
                        theme(panel.background = element_blank(), 
                              panel.border = element_blank(),
                              axis.line = element_line(color = "black"),
                              axis.ticks = element_line(color = "black"),
                              axis.ticks.length = unit(0.4, "lines"))
    } 
    else if(selected_method == "umap"){
      log_mat <- log2(raw_mat + 1) # log2(count+1)
      log_mat <- log_mat[apply(log_mat, 1, function(x) var(x) > 0), ]
      umap_res <- umap(t(log_mat), n_components=2, metric="euclidean")
      U_plot <- data.frame(UMAP1=as.numeric(umap_res[,1]),
                           UMAP2=as.numeric(umap_res[,2]),
                           Sample=rownames(umap_res))
      U_plot$group <- meta[[given_group]][match(U_plot$Sample, meta$Accession_ID)]
      d_plot <- ggplot(U_plot, aes(x = UMAP1, y = UMAP2, color = group, shape = group)) +
                        geom_point(size = 3) + # different shape
                        geom_text_repel(aes(label = Sample), size = 2, show.legend = FALSE) + # label do not overlap
                        theme_minimal() +
                        labs(x = "UMAP1", y = "UMAP2") +
                        theme(panel.background = element_blank(), 
                              panel.border = element_blank(),
                              axis.line = element_line(color = "black"),
                              axis.ticks = element_line(color = "black"),
                              axis.ticks.length = unit(0.4, "lines"))
    }
    
    values_basic$pca_plot <- d_plot
    # output
    output$distr_plot <- renderPlot({
      if(is.null(d_plot)) {
        return(NULL)
      }
      d_plot
    })
    
    resetLoadingButton("distr_button")
  }) # end of distr_button
  
  # Export the basic output
  output$export_basic <- downloadHandler(
    filename = function() {
      paste0("basic_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 7, height = 5)  # open PDF
      print(values_basic$depth_plot)
      print(values_basic$pca_plot)
      dev.off()  # close PDF
    },
    contentType = "application/pdf"
  )
  
  
  # Event 3 DE
  observeEvent(input$run_DE, {
    
    output$DE_finish <- renderUI({
      HTML("")
    })
    selected_method <- input$DE_method
    values_DE$method <- selected_method
    meta <- values_upload$meta
    raw_mat <- values_upload$mat
    log2FC <- input$log2FC
    q_val <- input$q_val
    design_form <- as.formula(paste("~", paste(input$DE_group, collapse=" + ")))
    
    if(selected_method == "deseq"){
      ddsMat <- DESeqDataSetFromMatrix(countData = raw_mat,
                                       colData = meta,
                                       design = design_form)
      
      keep <- rowSums(counts(ddsMat)) >= 10 # at least 10 reads total
      ddsMat <- ddsMat[keep,]
      # run
      dds <- DESeq(ddsMat) # including normalize
      values_DE$fit <- dds
      if(!is.null(dds)){
        output$DE_finish <- renderUI({
          HTML("<b>Differential analysis is complete!</b>")
        })
      }
    }
    
    get_name <- resultsNames(dds)
    updateSelectInput(session, "DE_name", choices=get_name)
    resetLoadingButton("run_DE")
  }) # end of run_DE
  
  
  # Event 4 Enriched genes
  observeEvent(input$vol_button, {
    
    # Debug for missing DE analysis
    if(is.null(values_DE$fit)){
      showNotification("Cannot find the Differential analysis results. Please run Differential analysis first.", type = "error")
      resetLoadingButton("vol_button")
      req(FALSE)
    }
    
    dds <- values_DE$fit
    de_name <- input$DE_name
    log2fc <- input$log2FC; q_val <- input$q_val
    values_DE$log2FC <- log2fc; values_DE$q_val <- q_val
    selected_org <- input$org_db_DE
    values_DE$org <- selected_org # reactive value
    org_map <- list(
      "Human" = "hsapiens_gene_ensembl",
      "Mouse" = "mmusculus_gene_ensembl",
      "Rat"   = "rnorvegicus_gene_ensembl"
    )
    org_db <- org_map[[selected_org]]
    
    if(values_DE$method == "deseq"){
      res <- results(dds, name=de_name)
      values_DE$res <- res
      MART <- useMart("ensembl", dataset = org_db)
      ensembl_ids <- rownames(res)
      genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = ensembl_ids,
                     mart = MART)
      
      symbol_map <- setNames(genes$hgnc_symbol, genes$ensembl_gene_id)
      gene_symbols <- symbol_map[rownames(res)]
      gene_symbols[is.na(gene_symbols)] <- rownames(res)[is.na(gene_symbols)]
      # ggplot volcano plot
      res$symbol <- gene_symbols
      vol_df <- data.frame(gene=rownames(res), log2FC=res$log2FoldChange, padj=res$padj, symbol=res$symbol)
      vol_df$Significant <- "No Significant"
      vol_df$Significant[vol_df$log2FC > log2fc & vol_df$padj < q_val] <- "Up"
      vol_df$Significant[vol_df$log2FC < -(log2fc) & vol_df$padj < q_val] <- "Down"
      # abs(log2FC) > log2fc + padj < q_val
      top_genes <- vol_df$gene[which(abs(vol_df$log2FC) > log2fc & vol_df$padj < q_val)]
      
      vol_df$tooltip <- paste0(
        "Gene: ", vol_df$gene, "\n",
        "Symbol: ", vol_df$symbol, "\n",
        "log2FC: ", signif(vol_df$log2FC, 3), "\n",
        "p_val: ", signif(vol_df$padj, 3)
      )
      
      vol_plot <- ggplot(vol_df, aes(log2FC, -log10(padj), text = tooltip))+
                          geom_point(aes(col=Significant))+
                          scale_color_manual(values=c("#0072B5", "grey", "#BC3C28"))+
                          labs(title = " ")+
                          geom_vline(xintercept=c(-(log2fc), log2fc), colour="black", linetype="dashed")+
                          geom_hline(yintercept = -log10(q_val),colour="black", linetype="dashed")+
                          theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
                          labs(x="log2(Fold-Change)",y="-log10(Pvalue)")+
                          theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
                          theme_bw()
      
      vol_plot_for_out <- ggplot(vol_df, aes(log2FC, -log10(padj)))+
                          geom_point(aes(col=Significant))+
                          scale_color_manual(values=c("#0072B5", "grey", "#BC3C28"))+
                          labs(title = " ")+
                          geom_vline(xintercept=c(-(log2fc), log2fc), colour="black", linetype="dashed")+
                          geom_hline(yintercept = -log10(q_val),colour="black", linetype="dashed")+
                          theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
                          labs(x="log2(Fold-Change)",y="-log10(Pvalue)")+
                          theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
                          theme_bw() +
                          geom_text_repel(
                            data = subset(vol_df, gene %in% top_genes),
                            aes(label = symbol),
                            size = 3)
                          
      }
      
      # store in reactive values
      values_DE$vol_plot <- vol_plot_for_out
      values_DE$DEG_table <- subset(vol_df, select = -tooltip)
      # output
      output$vol_plot <- renderPlotly({
        if(is.null(vol_plot)){
          return(NULL)
        }
        if (is.null(vol_plot)) return(NULL)
        ggplotly(vol_plot, tooltip = "text")
      })
      
      output$top_genes_table_ui <- renderUI({
        tagList(
          HTML("<h4><b>Enriched gene list</b></h4>"),
          div(DT::dataTableOutput("top_genes_table"))
        )
      })
      
      output$top_genes_table <- DT::renderDT({
          render_genes <- subset(vol_df, abs(log2FC) > log2fc & padj < q_val) %>%
                            mutate(across(where(is.numeric),
                                          ~ case_when(
                                            floor(.) == . ~ as.character(.),
                                            abs(.) > 1 ~ sprintf("%.2f", .),
                                            abs(.) < 0.001 ~ formatC(., format = "e", digits = 2),
                                            TRUE ~ sprintf("%.2f", .)
                                          )))
        render_genes <- subset(render_genes, select = -tooltip)
        datatable(render_genes, options=list(scrollX=T), rownames=FALSE)
      })
    
      resetLoadingButton("vol_button")
  }) # end of vol_button
  
  # control for vol_plot export
  observe({
    if (!is.null(values_DE$vol_plot)) shinyjs::show("export_volplot") else shinyjs::hide("export_volplot")
  })
  
  # Export the volcano plot output
  # output$export_volplot <- downloadHandler(
  #   filename = function() {
  #     paste0("volcano_", Sys.Date(), ".pdf")
  #   },
  #   content = function(file) {
  #     pdf(file, width = 7, height = 5)  # open PDF
  #     print(values_DE$vol_plot)
  #     dev.off()  # close PDF
  #   },
  #   contentType = "application/pdf"
  # )
  
  output$export_volplot <- downloadHandler(
    filename = function() {
      paste0("DEGs_", Sys.Date(), ".zip")
    },
    content = function(file) {
      tmpdir <- tempdir()
      temp <- file.path(tmpdir, as.integer(Sys.time())) # new folder under tmp
      dir.create(temp)
      f1 <- file.path(temp, "vol_plot.png")
      f2 <- file.path(temp, "DEGs.csv")
      
      ggsave(f1, plot=values_DE$vol_plot, width=6, height=4, dpi=300)
      write.csv(values_DE$DEG_table, file=f2, row.names=F)
      zip::zipr(zipfile=file, files=c("vol_plot.png", "DEGs.csv"), root=temp)
    },
    contentType = "application/zip"
  )
  
  # Event 5 pathway analysis
  output$dynamic_path_ui <- renderUI({
    if(input$path_db == 'GO'){
      tagList(
        fluidRow(
          column(6, numericInput("n_perm_go", label="nPerm", value=1000, min=10, max=9999)),
          column(6, numericInput("p_val_go", label="p_val cutoff", value=0.05, min=0, max=10))
        ),
        fluidRow(
          column(6, numericInput("minGGsize_go", label="minGGsize", value=3, min=1, max=9999)),
          column(6, numericInput("maxGGsize_go", label="maxGGsize", value=800, min=1, max=9999))
        ),
        fluidRow(
          column(6, selectInput("ont_go", label="Ontology", choices=c("BP", "MF", "CC", "ALL"))),
          column(6, numericInput("top_num", label="Show top genes (n=)", value=10, min=5, max=500))
        )
      )
    } else if(input$path_db == 'KEGG') {
      tagList(
        fluidRow(
          column(6, numericInput("n_perm_kk", label="nPerm", value=1000, min=10, max=9999)),
          column(6, numericInput("p_val_kk", label="p_val cutoff", value=0.05, min=0, max=10))
        ),
        fluidRow(
          column(6, numericInput("minGGsize_kk", label="minGGsize", value=3, min=1, max=9999)),
          column(6, numericInput("maxGGsize_kk", label="maxGGsize", value=800, min=1, max=9999))
        ),
        fluidRow(
          column(6, numericInput("top_num", label="Show top genes (n=)", value=10, min=5, max=500))
        )
      )
    }
  })
  
  observeEvent(input$path_button, {
    
    # Debug for missing DE analysis
    if(is.null(values_DE$res)){
      showNotification("Please run the Differential analysis and Enriched gene module first.", type="error")
      req(FALSE)
      resetLoadingButton("path_button")
    }
    
    DB <- input$path_db
    show_num <- input$top_num
    res <- values_DE$res
    selected_org <- values_DE$org
    
    # Debug for mismatch org_db
    if(input$org_db_DE != selected_org) {
      showNotification("Selected org database for DE was different from selected database in the pathway analysis.", type="error")
      req(FALSE)
      resetLoadingButton("path_button")
    }
    
    org_map <- list(
      "Human" = "org.Hs.eg.db",
      "Mouse" = "org.Mm.eg.db",
      "Rat"   = "org.Rn.eg.db"
    )
    org_db <- org_map[[selected_org]]
    tryCatch({
      suppressPackageStartupMessages(library(org_db, character.only = TRUE))
    }, error = function(e) {
      showNotification(paste0("Failed to load ", org_db, ": ", e$message,
                             "Please install the library ", org_db, "."), type = "error")
      req(FALSE)
      resetLoadingButton("path_button")
    })
    
    # initialize
    doseplot <- NULL
    path_res <- NULL
    gene_list <- NULL
    convert <- NULL
    gse <- NULL
    kk <- NULL
    
    # GSEA (signed fold change * -log10pvalue)
    res$fcsign <- sign(res$log2FoldChange)
    res$logP <- -log10(res$pvalue)
    res$metric <- res$fcsign * res$logP
    id_vec <- gsub("\\..*","", rownames(res))
    # log2 fold change
    original_gene_list <- res$metric
    # name the vector
    names(original_gene_list) <- id_vec
    # omit any NA values
    gene_list <- na.omit(original_gene_list)
    # sort the list in decreasing order (required for clusterProfiler)
    gene_list <- sort(gene_list, decreasing = TRUE)
    values_path$gene_list_ensembl <- gene_list
    values_path$gene_list_entrez <- gene_list
    
    convert <- data.frame(ENSEMBL = names(gene_list))
    get_db <- get(org_db)
    convert$ENTREZID <- mapIds(get_db,
                                keys = convert$ENSEMBL,
                                column = "ENTREZID",
                                keytype = "ENSEMBL",
                                multiVals = "first")
    
    if(DB == "GO") {
      n_perm <- input$n_perm_go
      p_val <- input$p_val_go
      minGGsize <- input$minGGsize_go
      maxGGsize <- input$maxGGsize_go
      ont <- input$ont_go
      # main
      go_gene_list <- values_path$gene_list_ensembl
      gse <- gseGO(geneList=go_gene_list, 
                   ont = ont, 
                   keyType = "ENSEMBL", 
                   nPerm = n_perm, 
                   minGSSize = minGGsize, 
                   maxGSSize = maxGGsize, 
                   pvalueCutoff = p_val, 
                   verbose = TRUE, 
                   OrgDb = org_db,
                   pAdjustMethod = "BH")
      # Dose
      if(nrow(gse@result) == 0){
        doseplot <- NULL; path_res <- NULL
        showNotification("Cannot find any enriched pathways.", type="warning")
      } else {
        doseplot <- dotplot(gse, showCategory=show_num, split=".sign", font.size=10, title="Gene Ontology") + facet_grid(.~.sign)
        path_res <- gse@result
        path_res <- path_res[, c("Description", "enrichmentScore", "qvalue")]
      }

      
    } else if(DB == "KEGG") {
      n_perm <- input$n_perm_kk
      p_val <- input$p_val_kk
      minGGsize <- input$minGGsize_kk
      maxGGsize <- input$maxGGsize_kk
      org_str_map <- list(
        "Human" = "hsa",
        "Mouse" = "mmu",
        "Rat"   = "rno"
      )
      org_str <- org_str_map[[selected_org]]
      # main
      kk_gene_list <- values_path$gene_list_entrez
      names(kk_gene_list) <- convert$ENTREZID
      kk <- gseKEGG(kk_gene_list, organism = org_str, 
                    keyType = "kegg", exponent = 1, 
                    nPerm = n_perm, minGSSize = minGGsize, 
                    maxGSSize = maxGGsize, pvalueCutoff = p_val, 
                    pAdjustMethod = "BH", verbose = TRUE, 
                    use_internal_data = FALSE, seed = FALSE)
      # Dose
      if(nrow(kk@result) == 0) {
        doseplot <- NULL; path_res <- NULL
        showNotification("Cannot find any enriched pathways.", type="warning")
      } else {
        doseplot <- dotplot(kk, showCategory=show_num, title="KEGG", split=".sign", font.size=10) + facet_grid(.~.sign)
        path_res <- kk@result
        path_res <- path_res[, c("Description", "enrichmentScore", "qvalue")]
      }

    }
    
    # store in reactive values
    values_path$path_dot <- doseplot
    values_path$path_table <- path_res
    
    # output
    output$top_path_table_ui <- renderUI({
      tagList(
        HTML("<h4><b>Enriched pathway list</b></h4>"),
        div(DT::dataTableOutput("top_path_table"))
      )
    })
    
    output$path_plot <- renderPlot({
      if(is.null(doseplot)) { return(NULL) }
      doseplot
    })
    
    output$top_path_table <- DT::renderDT({
      if(is.null(path_res)){ return(NULL) }
      render_res <- path_res %>%
        mutate(across(where(is.numeric),
                      ~ case_when(
                        floor(.) == . ~ as.character(.),
                        abs(.) > 1 ~ sprintf("%.2f", .),
                        abs(.) < 0.001 ~ formatC(., format = "e", digits = 2),
                        TRUE ~ sprintf("%.2f", .)
                      )))
      datatable(render_res, options=list(scrollX=T))
    })
    
    
    resetLoadingButton("path_button")
  }) # end of path_button
  
  # control for pathway export
  observe({
    if (!is.null(values_path$path_dot)) shinyjs::show("export_pathway") else shinyjs::hide("export_pathway")
  })
  
  # # Export the pathway plot output
  # output$export_pathway <- downloadHandler(
  #   filename = function() {
  #     paste0("pahtways_", Sys.Date(), ".pdf")
  #   },
  #   content = function(file) {
  #     pdf(file, width = 7, height = 5)  # open PDF
  #     print(values_path$path_dot)
  #     dev.off()  # close PDF
  #   },
  #   contentType = "application/pdf"
  # )
  
  output$export_pathway <- downloadHandler(
    filename = function() {
      paste0("enriched_pathways_", Sys.Date(), ".zip")
    },
    content = function(file) {
      tmpdir <- tempdir()
      temp <- file.path(tmpdir, as.integer(Sys.time()))
      dir.create(temp)
      f1 <- file.path(temp, "dotplot.png")
      f2 <- file.path(temp, "pathways.csv")
      
      ggsave(f1, plot=values_path$path_dot, width=8, height=10, dpi=300)
      write.csv(values_path$path_table, file=f2, row.names=F)
      zip::zipr(zipfile=file, files=c("dotplot.png", "pathways.csv"), root=temp)
    },
    contentType = "application/zip"
  )
  
} # server


# Run the application
shinyApp(ui=ui, server=server)
