library(shiny)
library(ggplot2)
library(plotly)
library(data.table)

# data preparation
load("data.rda")


# make long list format to accomodate both
plotPeptideComparison <- function(pepTraces_I, pepTraces_M, replicate = 1, uniprot_id = "protein_id"){
  
  pepTrace_ids_i <- pepTraces_I$trace_annotation[protein_id == uniprot_id, unique(id)]
  pepTraces_Is <- subset(pepTraces_I, trace_ids = pepTrace_ids_i)
  pepTraces_Is_long <- melt(pepTraces_Is$traces)
  setnames(pepTraces_Is_long, c("variable", "value"), c("fraction", "intensity"))
  pepTraces_Is_long$condition <- "Interphase"
  
  pepTrace_ids_m <- pepTraces_M$trace_annotation[protein_id == uniprot_id, unique(id)]
  pepTraces_Ms <- subset(pepTraces_M, trace_ids = pepTrace_ids_m)
  pepTraces_Ms_long <- melt(pepTraces_Ms$traces)
  setnames(pepTraces_Ms_long, c("variable", "value"), c("fraction", "intensity"))
  pepTraces_Ms_long$condition <- "Mitosis"
  
  pepTraces_IsMs_long <- rbind(pepTraces_Is_long, pepTraces_Ms_long)
  pepTraces_IsMs_long <- merge(pepTraces_IsMs_long, pepTraces_M$trace_annotation[, .(id, protein_id, Entry_name, Protein_names, Gene_names, protein_mw)],
                               by = "id", all.x = TRUE)
  pepTraces_IsMs_long
}

annotations <- names(trace_annotation_cum)
for (i in seq_along(annotations)){
  assign(annotations[i], trace_annotation_cum[[i]])
}

# server definition

shinyServer(function(input, output) {
  
  ## Generate Reactive Filter Value Field for UI, depending on filter column chosen
  output$fcolumnvalues <- renderUI({
    values <- sort(unique(get(input$fcolumn)))
    # values <- values[nchar(values)>0]
    selectizeInput("fvalue", "Search or select protein of interest", values, multiple = FALSE, options = list(maxOptions = 60000), selected = values[4533])
  })
  
  ## generate selected protein SEC traces plot
  output$protplot <- renderPlotly({
    
    # collect chromatograms for selection
    ######################################
    ######################################
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$fvalue),
                                      unique(protein_id)]
    # interphase
    target_id_traces_int_r1 <- toLongFormat(subset(prot_int_r1, trace_ids = target_id)$traces)
    target_id_traces_int_r1$replicate <- 1
    target_id_traces_int_r1$condition <- "Interphase"
    
    target_id_traces_int_r2 <- toLongFormat(subset(prot_int_r2, trace_ids = target_id)$traces)
    target_id_traces_int_r2$replicate <- 2
    target_id_traces_int_r2$condition <- "Interphase"
    
    target_id_traces_int_r3 <- toLongFormat(subset(prot_int_r3, trace_ids = target_id)$traces)
    target_id_traces_int_r3$replicate <- 3
    target_id_traces_int_r3$condition <- "Interphase"
    
    # mitosis
    target_id_traces_mit_r1 <- toLongFormat(subset(prot_mit_r1, trace_ids = target_id)$traces)
    target_id_traces_mit_r1$replicate <- 1
    target_id_traces_mit_r1$condition <- "Mitosis"
    
    target_id_traces_mit_r2 <- toLongFormat(subset(prot_mit_r2, trace_ids = target_id)$traces)
    target_id_traces_mit_r2$replicate <- 2
    target_id_traces_mit_r2$condition <- "Mitosis"
    
    target_id_traces_mit_r3 <- toLongFormat(subset(prot_mit_r3, trace_ids = target_id)$traces)
    target_id_traces_mit_r3$replicate <- 3
    target_id_traces_mit_r3$condition <- "Mitosis"
    
    target_id_traces <- rbind(target_id_traces_int_r1,
                              target_id_traces_int_r2,
                              target_id_traces_int_r3,
                              target_id_traces_mit_r1,
                              target_id_traces_mit_r2,
                              target_id_traces_mit_r3)
    
    target_id_traces <- merge(target_id_traces, trace_annotation_cum,
                              by.x = "id", by.y = "protein_id", all.y = FALSE, all.x = TRUE)
    
    target_id_traces[, monomer_fraction:=calibration_functions$MWtoSECfraction(protein_mw)]
    
    # PLOT
    ########################
    ########################
    lx.frc <- seq(5,(ncol(prot_int_r1$traces)-1),5)
    lx <- paste(lx.frc , round(calibration_functions$SECfractionToMW(lx.frc), 1) , sep = '(' )
    lx <- paste(lx, "", sep = ')' )
    
    p <- ggplot(target_id_traces, aes(x=fraction, y=intensity)) +
      ggtitle(paste(trace_annotation_cum[protein_id %in% target_id, .(protein_id, Entry_name, Gene_names)], collapse = " ")) 
      
    p <- p + geom_line(aes(group=paste(Gene_names, condition, replicate), color = condition), size = 2, alpha = 0.6) +
      scale_x_continuous(breaks = lx.frc , labels = lx) + xlab("SEC fraction number (apparent MW [kDa])") +
      ylab("Protein level SWATH-MS intensity (top2 peptide sum)") +
      scale_color_manual(values = c("black", "red")) 
    
    p <- p + theme_bw()
    
    if(input$show_monomers){
      p <- p + geom_point(aes(x = monomer_fraction, y = 0), shape = 23, color = "black", fill = "white", size = 3) 
    }
    
    if (input$logscale){
      p <- p + scale_y_log10()
    }
    
    if (input$remove_legend){
      p <- p +  theme(legend.position = "none")
    }
    
    ggplotly(p)
    dev.off()
    p
  })
  
  # peptide plot
  ## generate selected protein SEC traces plot
  output$pepplot <- renderPlotly({
    
    # collect chromatograms for selection
    ######################################
    ######################################
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$fvalue),
                                    unique(protein_id)]
    # PLOT
    ########################
    if(input$pepplot_replicate == 1){
      pepTraces_IsMs_long <- plotPeptideComparison(pep_int_r1, pep_mit_r1, replicate = 1, uniprot_id = target_id)
    }
    if(input$pepplot_replicate == 2){
      pepTraces_IsMs_long <- plotPeptideComparison(pep_int_r2, pep_mit_r2, replicate = 2, uniprot_id = target_id)
    }
    if(input$pepplot_replicate == 3){
      pepTraces_IsMs_long <- plotPeptideComparison(pep_int_r3, pep_mit_r3, replicate = 3, uniprot_id = target_id)
    }
    # plot
    lx.frc <- seq(5,(ncol(prot_int_r1$traces)-1),5)
    lx <- paste(lx.frc , round(calibration_functions$SECfractionToMW(lx.frc), 1) , sep = '(' )
    lx <- paste(lx, "", sep = ')' )
    
    p <- ggplot(pepTraces_IsMs_long, aes(fraction, intensity, group = id, color = id)) + geom_line(size = 1, alpha = 0.8) +
      facet_wrap(~condition, ncol = 1) + theme_bw() + ggtitle(unique(paste(pepTraces_IsMs_long$protein_id,
                                                                           pepTraces_IsMs_long$Entry_name,
                                                                           pepTraces_IsMs_long$Gene_names))) +
      scale_x_discrete(breaks = lx.frc, labels = lx) + xlab("SEC fraction number(apparent MW[kDa])")
    
    if (input$remove_legend){
      p <- p +  theme(legend.position = "none")
    }
    
    ggplotly(p)
    dev.off()
    p
  })
  
  output$table <- renderDataTable({
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$fvalue),
                                      unique(protein_id)]
    unique(merge(trace_annotation_cum[protein_id == target_id], up, by = "Entry_name", all.x = TRUE))
  })
  
  # Shift map plot(ly)s
  DiffCorrTableAvg <- as.data.table(as.data.frame(DiffCorrTable))
  DiffCorrTableAvg <-  unique(DiffCorrTableAvg[, .(protein_id, pearson_cor_avg, pearson_cor_avg_rank, difference_normalized_avg, difference_normalized_avg_rank, intensity_sum_log10_avg)])
  DiffCorrTableAvg <- merge(DiffCorrTableAvg, up, by.x = "protein_id", by.y = "Entry", all.x = TRUE)
  
  output$shiftmap_combined <- renderPlotly({
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$fvalue),
                                      unique(protein_id)]
    DiffCorrTableAvg$target <- DiffCorrTableAvg[[input$fcolumn]] %in% input$fvalue
    
    
    plot_ly(DiffCorrTableAvg, x = ~pearson_cor_avg, y = ~difference_normalized_avg,
                 type = 'scatter', mode = 'markers',  hoverinfo = 'text', alpha = 0.5,
            color = ~target, colors = c("black", "red"),  symbol = ~target, symbols = c("o", "x"),
                 text = ~paste('Uniprot_id: ', protein_id, 
                               '</br> Entry_name: ', Entry_name,
                               '</br> Gene_names: ', Gene_names,
                               '</br> Pathway: ', Pathway)) %>%
      layout(title = "Chromatogram Shift Score Map",
             xaxis = list(autorange = "reversed",
                          title = "Average Pearson Correlation I vs M"),
             yaxis = list(title = "Average Fractional Difference I vs M"))
  })
  
  output$shift_scores <- renderDataTable({
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$fvalue),
                                      unique(protein_id)]
    
    DiffCorrTable[protein_id == target_id, .(protein_id, Entry_name, Gene_names, replicate, pearson_cor, pearson_cor_rank, difference_normalized, difference_normalized_rank, shift_rank_sum)]
  })
  
  output$shiftmap_pearson <- renderPlotly({
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$fvalue),
                                      unique(protein_id)]
    DiffCorrTableAvg$target <- DiffCorrTableAvg[[input$fcolumn]] %in% input$fvalue
    
    plot_ly(DiffCorrTableAvg, x = ~pearson_cor_avg, y = ~intensity_sum_log10_avg,
            type = 'scatter', mode = 'markers',  hoverinfo = 'text', alpha = 0.5,
            color = ~target, colors = c("black", "red"),  symbol = ~target, symbols = c("o", "x"),
            text = ~paste('Uniprot_id: ', protein_id, 
                          '</br> Entry_name: ', Entry_name,
                          '</br> Gene_names: ', Gene_names,
                          '</br> Pathway: ', Pathway)) %>%
      layout(title = "Chromatogram Shift Score Map",
             xaxis = list(autorange = "reversed",
                          title = "Average Pearson Correlation I vs M"),
             yaxis = list(title = "Average Abundance (log10 cumulative Intensity I + M)"))
  })
  
  output$shiftmap_fracdiff <- renderPlotly({
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$fvalue),
                                      unique(protein_id)]
    DiffCorrTableAvg$target <- DiffCorrTableAvg[[input$fcolumn]] %in% input$fvalue
    
    target_id <- trace_annotation_cum[which(trace_annotation_cum[[input$fcolumn]] == input$fvalue),
                                      unique(protein_id)]
    
    plot_ly(DiffCorrTableAvg, x = ~difference_normalized_avg, y = ~intensity_sum_log10_avg,
            type = 'scatter', mode = 'markers',  hoverinfo = 'text', alpha = 0.5,
            color = ~target, colors = c("black", "red"),  symbol = ~target, symbols = c("o", "x"),
            text = ~paste('Uniprot_id: ', protein_id, 
                          '</br> Entry_name: ', Entry_name,
                          '</br> Gene_names: ', Gene_names,
                          '</br> Pathway: ', Pathway)) %>%
      layout(title = "Chromatogram Shift Score Map",
             xaxis = list(autorange = "reversed",
                          title = "Average Fractional Difference I vs M"),
             yaxis = list(title = "Average Abundance (log10 cumulative Intensity I + M)"))
  })
  
})
