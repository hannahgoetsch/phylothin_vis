#Shiny App - Phylogenetische-Baum-Analyse - Boschra Al Omari


library(shiny)
library(ape)


# Beispielbaum

base_path <- normalizePath(getwd(), mustWork = TRUE)

example_tree <- file.path(base_path, "ListeriaMonocytogenes.nwk")
example_removed <- file.path(base_path, "removed_ids_ListeriaMonocytogenes.txt")
example_clades <- file.path(base_path, "clades_ListeriaMonocytogenes.csv")

phylothin_path <- file.path(base_path, "phylothin-main/phylothin.R")
input_dir <- file.path(base_path, "phylothin_input")

output_dir <- file.path(input_dir, "phylothinoutput")



# Hilfsfunktionen

read_removed_file <- function(path) {
  readLines(path, warn = FALSE)
}

read_clades_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
#  lines <- gsub('^"|"$', '', lines)
#  lines <- gsub('""', '"', lines)
  
  df <- read.csv(
    text = lines,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"samples" %in% colnames(df)) {
    colnames(df)[2] <- "samples"
  }
  
  if (!"clade" %in% colnames(df)) {
    colnames(df)[3] <- "clade"
  }
  
  df <- df[, c("samples", "clade")]
  df$clade <- as.character(df$clade)
  
  return(df)
}

read_newick_safe <- function(path) {
  tree_text <- paste(readLines(path, warn = FALSE), collapse = "")
  
  if (!grepl(";$", tree_text)) {
    tree_text <- paste0(tree_text, ";")
    writeLines(tree_text, path)
  }
  
  read.tree(path)
}


# UI

ui <- fluidPage(
  
  titlePanel("Phylogenetischer Baum"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h3("Dateien auswählen"),
      
      radioButtons(
        "mode",
        "Datenquelle:",
        choices = c("Beispielbaum", "Eigener Baum"),
        selected = "Beispielbaum"
      ),
      
      conditionalPanel(
        condition = "input.mode == 'Eigener Baum'",
        fileInput("tree_file", "Eigenen Baum hochladen (.nwk)")
      ),
      
      textInput(
        "search_genome",
        "Genom suchen:",
        placeholder = "z.B. GCF_013282665"
      ),
      
      actionButton("run", "Analyse starten"),
      
      br(), br(),
      
      h4("Downloads"),
      downloadButton("download_reduced", "Reduzierten Baum herunterladen"),
      br(), br(),
      downloadButton("download_removed", "Removed IDs herunterladen"),
      br(), br(),
      downloadButton("download_clades", "Clusters herunterladen")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Originalbaum", plotOutput("originalPlot", height = "1000px")),
        tabPanel("Markierte Tips", plotOutput("markedPlot", height = "1000px")),
        tabPanel("Reduzierter Baum", plotOutput("reducedPlot", height = "1000px")),
        tabPanel("Alle Clusters", plotOutput("cladePlot", height = "1000px")),
        tabPanel("Große Clusters", plotOutput("bigCladePlot", height = "1000px")),
        tabPanel("Statistik", tableOutput("statsTable")),
        tabPanel("Vergleich", tableOutput("compareTable")),
        tabPanel("Cluster-Größen", plotOutput("cladeSizePlot", height = "700px"))
      )
    )
  )
)



# SERVER

server <- function(input, output, session) {
  
  tree <- reactiveVal(NULL)
  removed <- reactiveVal(NULL)
  clades <- reactiveVal(NULL)
  reduced <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    
    withProgress(message = "Analyse läuft...", value = 0, {
      
      incProgress(0.05, detail = "Analyse wird gestartet...")
      

      # Beispielbaum

      if (input$mode == "Beispielbaum") {
        
        incProgress(0.25, detail = "Beispielbaum wird geladen...")
        
        tr <- read_newick_safe(example_tree)
        rem <- read_removed_file(example_removed)
        cl <- read_clades_file(example_clades)
        
        incProgress(0.35, detail = "Baumdaten werden verarbeitet...")
        
        tree(tr)
        removed(rem)
        clades(cl)
        reduced(drop.tip(tr, rem))
        
        incProgress(0.35, detail = "Plots werden vorbereitet...")
      }
      

      # Eigener Baum

      else {
        
        req(input$tree_file)
        
        incProgress(0.10, detail = "Eigener Baum wird hochgeladen...")
        
        dir.create(input_dir, showWarnings = FALSE)
        
        tree_path <- file.path(input_dir, input$tree_file$name)
        file.copy(input$tree_file$datapath, tree_path, overwrite = TRUE)
        
        incProgress(0.10, detail = "Baum wird geprüft...")
        
        tr <- read_newick_safe(tree_path)
        
        if (is.null(tr)) {
          showNotification(
            "Der hochgeladene Baum konnte nicht gelesen werden.",
            type = "error"
          )
          return()
        }
        
        tree(tr)
        
        incProgress(0.20, detail = "PhyloThin wird ausgeführt...")
        
        old_wd <- getwd()
        setwd(base_path)
        
        system(
          paste(
            "Rscript",
            shQuote(phylothin_path),
            shQuote(input_dir),
            shQuote(basename(tree_path)),
            "no_PATHd8"
          )
        )
        
        setwd(old_wd)
        
        incProgress(0.20, detail = "PhyloThin-Ausgaben werden gesucht...")
        
        Sys.sleep(2)
        
        print("Dateien im Output-Ordner:")
        print(list.files(output_dir, full.names = FALSE))
        

        # Removed IDs laden

        removed_file <- list.files(
          output_dir,
          pattern = "removed_ids.*\\.txt$",
          full.names = TRUE,
          ignore.case = TRUE
        )
        
        if (length(removed_file) > 0) {
          
          rem <- read_removed_file(removed_file[1])
          rem <- trimws(rem)
          removed(rem)
          
        } else {
          
          kept_file <- list.files(
            output_dir,
            pattern = "kept_ids.*\\.txt$",
            full.names = TRUE,
            ignore.case = TRUE
          )
          
          if (length(kept_file) > 0) {
            
            kept <- readLines(kept_file[1], warn = FALSE)
            kept <- trimws(kept)
            
            rem <- setdiff(tr$tip.label, kept)
            removed(rem)
            
          } else {
            
            showNotification(
              "Keine removed_ids oder kept_ids Datei gefunden.",
              type = "error"
            )
            
            removed(character())
          }
        }
        
        incProgress(0.10, detail = "Reduzierter Baum wird geladen...")
        

        # Reduzierten Baum laden

        reduced_file <- list.files(
          output_dir,
          pattern = "reduced.*\\.nwk$|pruned.*\\.nwk$",
          full.names = TRUE,
          ignore.case = TRUE
        )
        
        if (length(reduced_file) > 0) {
          
          reduced(read_newick_safe(reduced_file[1]))
          
        } else {
          
          showNotification(
            "Kein reduzierter Baum gefunden. Reduced Tree wird aus removed_ids berechnet.",
            type = "warning"
          )
          
          reduced(drop.tip(tr, removed()))
        }
        
        incProgress(0.10, detail = "Clusters werden geladen...")
        

        # Clades laden

        clades_file <- list.files(
          output_dir,
          pattern = "clades.*\\.csv$",
          full.names = TRUE,
          ignore.case = TRUE
        )
        
        if (length(clades_file) > 0) {
          
          clades(read_clades_file(clades_file[1]))
          
        } else {
          
          showNotification(
            "Keine Clusters-Datei gefunden. Cluster-Plots können nicht angezeigt werden.",
            type = "warning"
          )
          
          clades(NULL)
        }
      }
      
      incProgress(1, detail = "Analyse abgeschlossen.")
    })
  })
  
  

  # Originalbaum

  output$originalPlot <- renderPlot({
    req(tree())
    
    tr <- tree()
    
    plot(
      tr,
      cex = 0.3,
      main = "Originalbaum"
    )
    
    if (input$search_genome != "") {
      searched <- grep(input$search_genome, tr$tip.label, ignore.case = TRUE)
      
      if (length(searched) > 0) {
        tiplabels(
          pch = 19,
          col = "blue",
          tip = searched,
          cex = 0.8
        )
      }
    }
  })
  
  

  # Markierte Tips

  output$markedPlot <- renderPlot({
    req(tree(), removed())
    
    tr <- tree()
    rem <- removed()
    
    removed_tips <- which(tr$tip.label %in% rem)
    
    plot(
      tr,
      cex = 0.3,
      main = "Originalbaum mit markierten Genomen"
    )
    
    tiplabels(
      pch = 4,
      col = "red",
      tip = removed_tips,
      cex = 0.6
    )
    
    if (input$search_genome != "") {
      searched <- grep(input$search_genome, tr$tip.label, ignore.case = TRUE)
      
      if (length(searched) > 0) {
        tiplabels(
          pch = 19,
          col = "blue",
          tip = searched,
          cex = 0.8
        )
      }
    }
  })
  
  

  # Reduzierter Baum

  output$reducedPlot <- renderPlot({
    req(reduced())
    
    tr <- reduced()
    
    plot(
      tr,
      cex = 0.4,
      main = "Reduzierter Baum"
    )
    
    if (input$search_genome != "") {
      searched <- grep(input$search_genome, tr$tip.label, ignore.case = TRUE)
      
      if (length(searched) > 0) {
        tiplabels(
          pch = 19,
          col = "blue",
          tip = searched,
          cex = 0.8
        )
      }
    }
  })
  
  

  # Alle Clades

  output$cladePlot <- renderPlot({
    req(reduced(), clades())
    
    tr <- reduced()
    cl <- clades()
    
    clade_map <- cl$clade[
      match(tr$tip.label, cl$samples)
    ]
    
    clade_map[is.na(clade_map)] <- "NA"
    
    clade_factor <- as.factor(clade_map)
    
    legend_labels <- levels(clade_factor)
    
    cols <- rainbow(length(legend_labels))
    names(cols) <- legend_labels
    
    tip_cols <- cols[clade_factor]
    
    plot(
      tr,
      tip.color = tip_cols,
      cex = 0.5,
      main = "Baum mit allen Clusters"
    )
    
    legend(
      "topright",
      legend = legend_labels,
      col = cols,
      pch = 19,
      bty = "o",
      cex = 0.6
    )
    
    if (input$search_genome != "") {
      searched <- grep(input$search_genome, tr$tip.label, ignore.case = TRUE)
      
      if (length(searched) > 0) {
        tiplabels(
          pch = 4,
          col = "black",
          tip = searched,
          cex = 0.9
        )
      }
    }
  })
  
  

  # Große Clades

  output$bigCladePlot <- renderPlot({
    req(reduced(), clades())
    
    tr <- reduced()
    cl <- clades()
    
    clade_map <- cl$clade[
      match(tr$tip.label, cl$samples)
    ]
    
    clade_map[is.na(clade_map)] <- "rest"
    
    tab <- table(cl$clade)
    big_clades <- names(tab[tab > 10])
    
    clade_factor2 <- clade_map
    clade_factor2[!clade_factor2 %in% big_clades] <- "rest"
    
    clade_factor2 <- as.factor(clade_factor2)
    
    legend_labels2 <- levels(clade_factor2)
    
    cols2 <- rainbow(length(legend_labels2))
    names(cols2) <- legend_labels2
    
    cols2["rest"] <- "grey70"
    
    tip_cols2 <- cols2[clade_factor2]
    
    plot(
      tr,
      tip.color = tip_cols2,
      cex = 0.5,
      main = "Große Clusters"
    )
    
    keep <- legend_labels2 != "rest"
    
    legend(
      "topright",
      legend = c(legend_labels2[keep], "Die restlichen Clusters"),
      col = c(cols2[keep], "grey70"),
      pch = 19,
      bty = "n",
      cex = 0.7
    )
    
    if (input$search_genome != "") {
      searched <- grep(input$search_genome, tr$tip.label, ignore.case = TRUE)
      
      if (length(searched) > 0) {
        tiplabels(
          pch = 4,
          col = "black",
          tip = searched,
          cex = 0.9
        )
      }
    }
  })
  
  

  # Statistik-Tab

  output$statsTable <- renderTable({
    req(tree(), reduced(), removed())
    
    cl <- clades()
    
    original_n <- length(tree()$tip.label)
    reduced_n <- length(reduced()$tip.label)
    removed_n <- length(removed())
    reduction_percent <- round((removed_n / original_n) * 100, 2)
    
    if (!is.null(cl)) {
      clade_count <- length(unique(na.omit(cl$clade)))
      na_count <- sum(is.na(cl$clade) | cl$clade == "NA")
    } else {
      clade_count <- NA
      na_count <- NA
    }
    
    data.frame(
      Kennzahl = c(
        "Anzahl Genome im Originalbaum",
        "Anzahl entfernte Genome",
        "Anzahl Genome im reduzierten Baum",
        "Reduktionsrate in %",
        "Anzahl berechneter Clusters",
        "Genome ohne Cluster-Zuordnung"
      ),
      Wert = c(
        original_n,
        removed_n,
        reduced_n,
        reduction_percent,
        clade_count,
        na_count
      )
    )
  })
  
  

  # Vergleich Original vs. Reduziert

  output$compareTable <- renderTable({
    req(tree(), reduced(), removed())
    
    original_n <- length(tree()$tip.label)
    reduced_n <- length(reduced()$tip.label)
    removed_n <- length(removed())
    
    data.frame(
      Kategorie = c("Genome", "Entfernte Genome", "Verbleibende Genome"),
      Originalbaum = c(original_n, 0, original_n),
      Reduzierter_Baum = c(reduced_n, removed_n, reduced_n)
    )
  })
  
  

  # Clade-Größendiagramm

  output$cladeSizePlot <- renderPlot({
    req(clades())
    
    cl <- clades()
    
    clade_values <- cl$clade
    clade_values[is.na(clade_values)] <- "NA"
    
    tab <- sort(table(clade_values), decreasing = TRUE)
    
    barplot(
      tab,
      las = 2,
      cex.names = 0.7,
      main = "Größe der Clusters",
      xlab = "Cluster",
      ylab = "Anzahl Genome"
    )
  })
  
  

  # Downloads

  output$download_reduced <- downloadHandler(
    filename = function() {
      "reduced_tree.nwk"
    },
    content = function(file) {
      req(reduced())
      write.tree(reduced(), file = file)
    }
  )
  
  output$download_removed <- downloadHandler(
    filename = function() {
      "removed_ids.txt"
    },
    content = function(file) {
      req(removed())
      writeLines(removed(), con = file)
    }
  )
  
  output$download_clades <- downloadHandler(
    filename = function() {
      "clades.csv"
    },
    content = function(file) {
      req(clades())
      write.csv(clades(), file = file, row.names = FALSE)
    }
  )
}



# APP STARTEN

shinyApp(ui = ui, server = server)
