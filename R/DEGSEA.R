# Cargamos las librerias necesarias
library(dplyr)
library(ggplot2)
library(data.table)
library(MatrixGenerics)
library(matrixStats)
library(tximport)
library(DESeq2)
# library(btools)
library(WGCNA)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(phyloseq)
library(EnhancedVolcano)
library(ggpubr)
library(clusterProfiler)
library(ashr)
library(devtools)
library(ontologyIndex)

##################################
#           FUNCIONES            #
##################################

#' @export
Notes2GO <- function(Notes, godb, column_GO, column_ID, separator, level=c("Gene","Tx")){

  separator <- paste0("\\",separator)

  if (level == "Gene") {
    GO_df <- Notes %>%
      dplyr::select(.data[[column_GO]], .data[[column_ID]]) %>%
      separate_rows(.data[[column_GO]], sep = separator) %>%
      mutate(gene = gsub("\\..*", "", .data[[column_ID]])) %>% # Elimina todo lo que vaya detras de un punto
      select(.data[[column_GO]], gene) %>%
      distinct() %>%
      drop_na()
    colnames(GO_df) <- c("term", "gene")
  }else{
    GO_df <- Notes %>%
      dplyr::select(.data[[column_GO]], .data[[column_ID]]) %>%
      separate_rows(.data[[column_GO]], sep = separator) %>%
      mutate(gene = .data[[column_ID]]) %>%
      select(.data[[column_GO]], gene) %>%
      distinct() %>%
      drop_na()
    colnames(GO_df) <- c("term", "gene")
  }

  # Se descarga de manera manual la DB de go.obo: https://geneontology.org/docs/download-ontology/

  # Se saca la ontologia y se crean tablas de interconversion para las funciones de clusterProfiler
  ontology <- get_ontology(file = godb,
                           propagate_relationships = "is_a",
                           extract_tags = "everything",
                           merge_equivalent_terms = TRUE)


  GO_df_term <- GO_df %>%
    mutate(name = ontology$name[term]) %>%
    select(c(term, name)) %>%
    distinct() %>%
    drop_na() %>%
    filter(!grepl("obsolete", name))

  GO_df <- GO_df %>%
    filter(term %in% GO_df_term$term)

  write.table(GO_df, file = "ter2geneGO.tsv", quote = FALSE, sep = "\t")
  write.table(GO_df_term, file = "ter2nameGO.tsv", quote = FALSE, sep = "\t")

  return(list(term2geneGO = GO_df, term2nameGO = GO_df_term))

}

#' @export
Deseq2_Vulcano <- function(deseq2_table, title_plot, plot_outname, subtitle){

  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  report = deseq2_table
  # Formateo preventivo (asumimos_guion_bajo_y_versus)
  title_plot = str_replace(title_plot, "DESeq2_", "")
  title_plot = str_replace(title_plot, "_vs_", " vs ")
  title_plot = str_replace(title_plot, "_plot", "")

  plot = EnhancedVolcano(report,
                         title = title_plot,
                         lab = rownames(report),
                         selectLab = FALSE,
                         titleLabSize = 15,
                         subtitle = subtitle,
                         subtitleLabSize = 15,
                         captionLabSize = 15,
                         x = "log2FoldChange",
                         y = "padj",
                         FCcutoff = 1,
                         labSize = 7,
                         # legend = NA,
                         boxedLabels = TRUE,
                         pointSize = 2.0,
                         # axisLabSize = 10,
                         max.overlaps = 5,
                         legendPosition = "none",
                         gridlines.major = TRUE,
                         gridlines.minor = TRUE,
                         drawConnectors = TRUE,
                         widthConnectors = 0.5,
                         pCutoff = 0.05)

  ggsave(plot = plot, filename = plot_outname, height = 5, width = 4, dpi = 300)

}


#' @export
Functional_Enrichment <- function(deseq2_table, out_name, term2gene, term2name, outdir, immune = FALSE){
  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  report = deseq2_table
  title = out_name
  title = str_replace(title, "DESeq2_", "")
  title = str_replace(title, "_plot", "")

  reportf = unique(report)

  ############################################################################
  # GSEA
  # This time I am not filtering the changes
  rownames(reportf) <- reportf$ID
  gsea_list <- reportf %>%
    # dplyr::filter(log2FoldChange >= 1 & padj <= 0.05) %>%
    dplyr::arrange(desc(log2FoldChange))

  gsea_input <- gsea_list %>%
    dplyr::select(log2FoldChange) %>%
    unlist() %>%
    as.vector()
  names(gsea_input) <- gsub("mRNA:","" ,rownames(gsea_list)) # De precaucion

  # Gene enrichment analysis using the whole dataset sorted by FC
  enrichment_gsea <- GSEA(geneList = gsea_input,
                          TERM2GENE = term2gene,
                          TERM2NAME = term2name,
                          minGSSize = 10,
                          maxGSSize = 500,
                          eps = 1e-10,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH")
  #save the enrichment result
  enr <- enrichment_gsea@result
  enr$core_enrichment <- gsub("/","|", enr$core_enrichment)
  enr$leading_edge <- NULL # No lo necesitamos y descuadra las columnas al abrar como csv
  enr$Description <- gsub(",",";", enr$Description)
  write.csv(file = paste0(outdir, "/", title, "_GSEA.csv"),
            x = enr, quote = FALSE, row.names = FALSE)

  if (any(enrichment_gsea@result$p.adjust <= 0.05)){
    p <- ridgeplot(enrichment_gsea,
                   core_enrichment= TRUE, # Cambiado, en FALSE no existe
                   fill="p.adjust",
                   orderBy = "NES",
                   showCategory=100) +
      ggtitle("Ridge plot for GSEA")

    ggsave(filename = paste0(outdir, "/", title, "_GSEA_ridgeplot.pdf"),
           plot =  p,  dpi = 300, width = 21, height = 70, units = "cm")
  }
  # Filtering by immune terms
  if (immune == TRUE) {
    enrichment_gsea_immuno <- `enrichment_gsea@result`
    enrichment_gsea_immuno <- enrichment_gsea_immuno[grepl("defens|immun| amp |antimicrobial|infection|toll|imd",
                                                           enrichment_gsea_immuno$Description,
                                                           ignore.case = TRUE),]

    enrichment_gsea_immuno <- enrichment_gsea_immuno[!grepl("t cell|b cell|lymphocyte|natural killer|mast cell|recombination",
                                                            enrichment_gsea_immuno$Description,
                                                            ignore.case = TRUE),]

    enrichment_gsea_immuno <- enrichment_gsea_immuno %>% filter(p.adjust < 0.05)

    if (nrow(enrichment_gsea_immuno) > 0) {
      write.csv(file = paste0(outdir, "/", title,"_Immuno_GSEA.csv"), x = enrichment_gsea_immuno)
    }
  }

}


#' @export
DE_pairwise <- function(mifactor, muestras, objeto_dds, subdir = FALSE, prefix = FALSE,
                        immune = FALSE, tipo=c("factor","subfactor"),
                        ter2nameGO, ter2geneGO, descriptions=FALSE){

  if (tipo == "factor") {
    # Creamos una carpetas de salida (si no es subfactor)
    outpath <- paste0(base::getwd(), "/", mifactor, "/")
    base::dir.create(outpath)
    outpath_vulcano <- paste0(outpath, "Vulcano_plots")
    base::dir.create(outpath_vulcano)
    outpath_tables <- paste0(outpath, "DESeq2_tables")
    base::dir.create(outpath_tables)
    outpath_GSEA <- paste0(outpath, "GSEA")
    base::dir.create(outpath_GSEA)
    if (immune == TRUE) {
      outpath_GSEA_Immuno <- paste0(outpath, "GSEA_Immuno")
      base::dir.create(outpath_GSEA_Immuno)
    }
  }else{
    # Creamos una carpetas de salida (si no es subfactor)
    outpath <- paste0(base::getwd(), "/", subdir, "/")
    base::dir.create(outpath)
    outpath_vulcano <- paste0(outpath, "Vulcano_plots")
    base::dir.create(outpath_vulcano)
    outpath_tables <- paste0(outpath, "DESeq2_tables")
    base::dir.create(outpath_tables)
    outpath_GSEA <- paste0(outpath, "GSEA")
    base::dir.create(outpath_GSEA)
    if (immune == TRUE) {
      outpath_GSEA_Immuno <- paste0(outpath, "GSEA_Immuno")
      base::dir.create(outpath_GSEA_Immuno)
    }
  }

  if (immune == TRUE) {
    # Inicializamos el dataframe para las listas de respuesta humoral
    Inmuno <- data.frame(matrix(ncol = 9, nrow = 0))
    x <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","Group","TxID")
    colnames(Inmuno) <- x
  }

  # Establecemos los grupos y las comparciones en funcion del factor
  lev <- levels(as.factor(muestras[,mifactor])) # get the variables
  L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

  for (j in 1:length(L.pairs)) {

    # Extraccion de los 2 factores de cada pareja posible
    factor1 = L.pairs[[j]][1]
    factor2 = L.pairs[[j]][2]

    # Creacion del nombre de la comparcion e inicalizacion del dataframe de resultados
    nam <- paste0("DESeq2_",factor1,"_vs_",factor2)
    if (tipo == "factor") {
      plot_n <- paste0(nam,"_plot")
    }else{
      plot_n <- paste0(nam, "_", prefix, "_plot")
    }

    dds_df <- data.frame(matrix(ncol = 7, nrow = 0))
    x <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","rank")
    colnames(dds_df) <- x

    # Aplicacion del test de Wald de DESeq2 con correccion "ash" del log2FoldChange
    factor1vs2 <- results(objeto_dds, contrast=c(mifactor, factor1, factor2), alpha = 0.05)
    factor1vs2 <- lfcShrink(objeto_dds, contrast=c(mifactor, factor1, factor2), res=factor1vs2, type = "ashr") #type = "normal") # Este tipo funciona
    factor1vs2 = data.frame(factor1vs2)
    dds_df = rbind(dds_df, factor1vs2)
    dds_df <- data.frame(dds_df)
    rownames(dds_df) <- gsub("mRNA:", "", rownames(dds_df)) # Eliminacion preventiva de cabecera rara htseq
    plot_n <- as.character(plot_n)

    # Filtrado por significancia
    reportf = unique(dds_df)
    reportf = reportf %>% filter(padj < 0.05)
    reportf = reportf %>% filter(abs(log2FoldChange) >= 0.6)
    subtitle_plot <- paste0("DEGs: ", nrow(reportf))

    # Test para ver si el nombre de los ids se aguanta
    reportf2 <- reportf
    reportf2$ID <- rownames(reportf)
    reportf2$RefGroup <- rep(factor2)
    reportf2$StuGroup <- rep(factor1)
    reportf2 <- reportf2[,c("ID","StuGroup","RefGroup","baseMean","log2FoldChange","padj")]
    # Usamos un try aqui
    try(reportf2 <- base::merge(reportf2, descriptions, by="ID"))

    if (immune == TRUE) {

      # Filtrado por listas IMD/Toll/AMPs

      selected_tx <- c()
      selected_tx <- c(selected_tx, rownames(reportf[grepl("_", rownames(reportf)),]))
      report_Inmuno <- reportf[grepl("_", rownames(reportf)),] # Usamos el patron "_"
      # report_Inmuno$Group <- rep(title_plot, nrow(report_Inmuno))

      # Usaremos este filtrado de lo que nos interesa para señalar los DEGs inmunos en el volcano plot
      report_Inmuno = unique(report_Inmuno)
      report_Inmuno$TxID <- rownames(report_Inmuno)
      Inmuno <- rbind(Inmuno, report_Inmuno)
    }


    # Creamos un Vulcano Plot
    Vulcname <- paste0(plot_n, ".png")
    Deseq2_Vulcano(deseq2_table = dds_df, title_plot = plot_n,
                   plot_outname = paste0(outpath_vulcano,"/",Vulcname), subtitle = subtitle_plot)

    # Creamos la tabla DESeq2 con resultados
    write.table(reportf2, file = paste0(outpath_tables,"/", gsub("_plot","_table", plot_n)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Analisis GSEA con la funcion definida
    Functional_Enrichment(deseq2_table = reportf2, out_name = plot_n,
                          term2name = ter2nameGO, term2gene = ter2geneGO,
                          outdir = outpath_GSEA)

  }

  # Creamos el archivo con los datos de las listas IMD/Toll/AMPs
  #write.table(Inmuno, file = paste0(outpath,"/Inmuno_DESeq2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}



Disect2 <- function(samples,
                    samplename = "sampleName",
                    filename = "fileName",
                    path = "path",
                    factor1, subfactor = FALSE,
                    feature=c("Gene","Tx"),
                    countool=c("kallisto","htseqcounts"),
                    term2gene, term2name,
                    directory,
                    descriptions_available=FALSE,
                    immune = FALSE){

  if (countool == "kallisto") {

    # Reordenamos el df considerando que al pasar los files en tximport, lo ordenara alfanumericamente
    rownames(samples) <- samples[[samplenames]]
    samples <- samples[order(samples[[path]], decreasing = FALSE),] # Nos aseguramos de mantener el orden

    # Pre-Importamos los conteos de h5
    txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE) # Evitamos meter tx2gene

    if (feature == "Gene") {
      # PARTE PARA AGREGAR POR GENES, NO SE HACE EN ESTE ANALISIS
      tx2gene <- data.frame(names(txi.kallisto$abundance[,1]), gsub("t[0-9].*","",(names(txi.kallisto$abundance[,1]))))
      colnames(tx2gene) <- c("txID","geneID")
      tx2gene$geneID <- gsub("_sm","_smt3", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3
      txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

    }

    for (f in factor1) {
      if (subfactor == FALSE) {
        # Creacion del objeto DESeq2
        dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = as.formula(paste0("~ ", f)))
        dds <- DESeq(dds)
        DE_pairwise(mifactor=f, muestras=samples, objeto_dds=dds, tipo = "factor",
                    ter2nameGO = term2name, ter2geneGO = ter2geneGO, descriptions = descriptions_available)
      }else{
        # Creamos los ficheros asqui directamente
        outpath <- paste0(base::getwd(), "/", f, "/")
        base::dir.create(outpath)
        outpath_vulcano <- paste0(outpath, "Vulcano_plots")
        base::dir.create(outpath_vulcano)
        outpath_tables <- paste0(outpath, "DESeq2_tables")
        base::dir.create(outpath_tables)
        outpath_GSEA <- paste0(outpath, "GSEA")
        base::dir.create(outpath_GSEA)
        if (immune == TRUE) {
          outpath_GSEA_Immuno <- paste0(outpath, "GSEA_Immuno")
          base::dir.create(outpath_GSEA_Immuno)
        }
        for (j in unique(samples[[f]])) {
          subsamples <- samples[grepl(j, samples[[f]]),]
          sf <- paste(subsamples$sample, collapse="|")
          subfiles <- files[grepl(sf, files)]
          sub_txi.kallisto <- tximport(subsamples$path, type = "kallisto", txOut = TRUE)
          des <- as.formula(paste0("~ ", subfactor))
          dds <- DESeqDataSetFromTximport(sub_txi.kallisto, colData = subsamples, design = des)
          dds <- DESeq(dds)
          DE_pairwise(mifactor = subfactor, muestras = subsamples, objeto_dds = dds,
                      tipo = "subfactor", subdir = f, prefix = j,
                      ter2nameGO = term2name, ter2geneGO = ter2geneGO, descriptions = descriptions_available)
        }
      }
    }
  }



  if (countool == "htseqcounts") {

    # Adding the appropiate order for htseqcount-deseq2 object creation
    sampleTable <- samples
    rownames(sampleTable) <- sampleTable$sampleName

    for (f in factor1) {
      if (subfactor == FALSE) {
        # sampleTable <- sampleTable[,c(samplename, filename, f)]

        # Hacemos un objeto DESeq2 con htseqcount y apartir de ahi reordenamos el assay
        des <- as.formula(paste0("~ ", f))
        dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                          directory = directory,
                                          design = des)
        dds <- DESeq(dds)
        DE_pairwise(mifactor=f, muestras=sampleTable, objeto_dds=dds, tipo = "factor",
                    ter2nameGO = term2name, ter2geneGO = ter2geneGO, descriptions = descriptions_available)

      }else{
        # Creamos los ficheros asqui directamente
        for (j in unique(sampleTable[[f]])) {
          subsamples <- sampleTable[grepl(j, sampleTable[[f]]),]
          subsamples[[f]] <- as.factor(subsamples[[f]])
          des <- as.formula(paste0("~ ", subfactor))
          dds <- DESeqDataSetFromHTSeqCount(sampleTable = subsamples,
                                            directory = directory,
                                            design = des)
          dds <- DESeq(dds)
          DE_pairwise(mifactor = subfactor, muestras = subsamples, objeto_dds = dds, immune = FALSE,
                      tipo = "subfactor", subdir = f, prefix = j,
                      ter2nameGO = term2name, ter2geneGO = ter2geneGO, descriptions = descriptions_available)
        }
      }
    }
  }
}



