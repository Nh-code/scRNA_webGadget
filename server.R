#Author : Nh_code
#Date: 2021-12-27
#Filename: scRNA_webGadget.R
#Version: 1.0
#Contact_author: zhuzhiyong71@gmail.com
#Description: Deployment of some single-cell processes under the framework of shiny
################################################################################################
library(shiny)
library(shinyWidgets)
library(Seurat)
library(tidyverse)
library(shinycustomloader)
library(DT)
library(ggplot2)
library(threejs)
library(rvest)
library(biomaRt)
library(echarts4r)

function(input, output, session) {

  options(shiny.maxRequestSize=100000*1024^2)
  warnings('off')
    
    upload_myfile <- eventReactive(input$my_file,{
        rowdata <- readRDS(file = input$my_file$datapath)
        return(rowdata)
    })
    
    scRNA <- eventReactive(c(input$submit,input$refresh),{
        all <- upload_myfile()
        if (input$myIQR) {
            Q1 <- quantile(all$nFeature_RNA,0.25)
            Q3 <- quantile(all$nFeature_RNA,0.75)
            IQR <- Q3 - Q1
            upper <- Q3+1.5*IQR
            lower <- Q1-1.5*IQR
        }else{
            upper <- input$upper
            lower <- input$lower
        }
        all <- subset(all, subset = nFeature_RNA > lower & nFeature_RNA < upper)
        all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
        all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
        all.genes <- unique(rownames(all))
        all <- ScaleData(all, features = all.genes)
        all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
        dim.Usage <- seq( input$dims[1],input$dims[2] )
        all <- FindNeighbors( all, dims = dim.Usage )
        all <- FindClusters( all, resolution = input$res )
        all <- RunUMAP( all, dims = dim.Usage )
        all <- RunTSNE( all,dims = dim.Usage )
        return(all)
    })
    ################################### PAGE - 1 #######################################################
    
    output$summary <- shiny::renderUI({
        input$submit
        myEC <- isolate(scRNA())
        sumNfeature <- summary(myEC$nFeature_RNA)
        sumNcounts <- summary(myEC$nCount_RNA)
        CellNum <- paste0("<p style=\"color:#FF4500;font-size: 150%;\">CellNum: ",dim(myEC)[2],"</p>")
        nFeatures_range <- paste0("<p style=\"font-size: 125%;font-weight: bold\">nFeatures_range:   [  ",sumNfeature[1]," ~ ",sumNfeature[6],"  ]      Mean = ",sumNfeature[4],"</p>")
        nCounts_range <- paste0("<p style=\"font-size: 125%;font-weight: bold\">nCounts_range:   [  ",sumNcounts[1]," ~ ",sumNcounts[6],"  ]      Mean = ",sumNcounts[4],"</p>")
        HTML(paste(CellNum,nFeatures_range,nCounts_range))

    })

    output$metadata <- DT::renderDT({
        input$submit
        myEC <- isolate(scRNA())
        DT::datatable(myEC@meta.data,
                      colnames = c("CellBarcode" = 1),
                      options = list(
                        searchHighlight= T,
                        autoWidth = T,
                        columnDefs = list(list(className='dt-center',width = '100px', targets = c(1,6))),
                        pageLength = 10,
                        scrollX = T,
                        scrollY = T,
                        initComplete = JS(
                          "function(settings, json) {",
                          "$(this.api().table().header()).css({'background-color': '#00aeef', 'color': '#fff'});",
                          "}")
                      )
                      
        )
    })

    output$plot_1 <- shiny::renderPlot({
      input$submit
      myEC <- isolate(scRNA())
      myEC@meta.data$orig.ident <- "myObject"
      VlnPlot(myEC, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident",ncol = 2,pt.size = 0) 

    })
    
    ################################### PAGE - 2 #######################################################
   
     output$plot_pca <- threejs::renderScatterplotThree({
      input$refresh
      myEC <- isolate(scRNA())
      my_PCA_compo <- as.data.frame(myEC@reductions$pca@cell.embeddings)
      N <- length(levels(myEC@meta.data$seurat_clusters))
      x.usage <- my_PCA_compo[[input$PCA_X]]
      y.usage <- my_PCA_compo[[input$PCA_Y]]
      z.usage <- my_PCA_compo[[input$PCA_Z]]
      axi_label <- c(input$PCA_X,input$PCA_Z,input$PCA_Y)
      threejs::scatterplot3js(x = x.usage,y = y.usage,z = z.usage,
                              x.ticklabs =NA ,
                              y.ticklabs =NA ,
                              z.ticklabs =NA ,
                              grid = T,
                              signif = 10,
                              cex.symbols = 0.3,
                              pch = "@" ,
                              stroke = NULL,
                              axisLabels = axi_label,
                              col = rainbow(N)[myEC@meta.data$seurat_clusters]
      )
    })
    
    output$UMAP <- shiny::renderPlot({
      input$refresh
      myEC <- isolate(scRNA())
      DimPlot(myEC,reduction = "umap",label = TRUE)
      
    },height = 350)

    output$tSNE <- shiny::renderPlot({
      input$refresh
      myEC <- isolate(scRNA())
      DimPlot(myEC,reduction = "tsne",label = TRUE)
    },height = 350)
    
    ################################### PAGE - 3 #######################################################
   
     modify_vlnplot <- function(obj, feature, pt.size  , plot.margin = unit(c(0, 1, 0, 1), "cm"),...) {
      p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
        xlab("") + ylab(feature) + ggtitle("") + 
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_line(),
              axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
              plot.margin = plot.margin ) 
      return(p)
    }
    #color theme
    my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                    '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                    '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                    '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                    '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                    '#968175')
    ## wrap_vln function
    StackedVlnPlot <- function(obj, features , pt.size , title , plot.margin = unit(c(-0.75, 3, -0.75, 3), "cm"), ...) {
      if (missing(title)) title <- ""
      plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x,pt.size,...))
      plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
        theme(axis.text.x=element_text(), axis.ticks.x = element_line())
      plot_list[[1]]<- plot_list[[1]] +ggtitle(title)+theme(plot.title = element_text(size = rel(1),hjust = 0.5,vjust = 1.5))
      p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
      return(p)
    }
    NCBI_geninfo <- function(symbol){
      feature <- symbol
      mart <- useMart("ensembl","hsapiens_gene_ensembl")
      entrze_ID <- getBM(attributes = c('entrezgene_id','external_gene_name'),
                         filters = "external_gene_name",values = feature,
                         mart =mart)
      NCBI_url <- paste("https://www.ncbi.nlm.nih.gov/gene/",entrze_ID$entrezgene_id,sep="")
      NCBI_webpage <- rvest::read_html(NCBI_url)
      
      summary <- NCBI_webpage %>% html_nodes(xpath='//*[@id="summaryDl"]/dd[10]/text()') %>% html_text()
      Also_known_as<- NCBI_webpage %>% html_nodes(xpath='//*[@id="summaryDl"]/dd[9]/text()') %>% html_text()
      Organism <- NCBI_webpage %>% html_nodes(xpath='//*[@id="summaryDl"]/dd[7]/a/text()')  %>% html_text()
      Official_Full_Name <- NCBI_webpage %>% html_nodes(xpath='//*[@id="summaryDl"]/dd[2]/text()') %>% html_text()
      Official_Symbol <- NCBI_webpage %>% html_nodes(xpath='//*[@id="summaryDl"]/dd[1]/text()') %>% html_text()
      
      my_Also_known_as <- paste("<span style=\"color:#515151;font-size: 150%;font-weight: bolder\">Also known as:&nbsp;&nbsp;</span><span style=\"color:#515151;font-size: 150%\">",Also_known_as,"</span>",sep = "")
      my_Organism <- paste("<span style=\"color:#515151;font-size: 150%;font-weight: bolder\">Organism :&nbsp;&nbsp;</span><span style=\"color:#515151;font-size: 150%\">",Organism,"</span>",sep = "")
      my_Official_Full_Name <- paste("<span style=\"color:#515151;font-size: 150%;font-weight: bolder\">Official Full Name:&nbsp;&nbsp;</span><span style=\"color:#515151;font-size: 150%\">",Official_Full_Name,"</span>",sep = "")
      my_Official_Symbol <- paste("<span style=\"color:#515151;font-size: 150%;font-weight: bolder\">Official Symbol:&nbsp;&nbsp;</span><span style=\"color:#515151;font-size: 150%\">",Official_Symbol,"</span>",sep = "")
      
      my_summary <- paste("<span style=\"color:#515151;font-size: 150%;font-weight: bolder\">Summary:&nbsp;&nbsp;</span><span style=\"color:#0073b6;font-size: 120%\"><i>",
                          summary,
                          "</i></span>",sep = "")
      text_1 <- paste(my_Official_Symbol,my_Official_Full_Name,my_Organism,my_Also_known_as,sep = "</br>")
      return(paste(text_1,my_summary,sep = "</br></br>"))
    }
    cancer_npmi <- function(symbol){
      feature <- symbol
      my_url <- paste0("http://chat.lionproject.net/chartdata?q=",toupper(feature),"&measure=npmi")
      my_npmi <- my_url %>% httr::GET(config = httr::config(ssl_verifypeer = FALSE)) %>%  read_html()  %>% html_text()
      my_npmi_vc <- unlist(strsplit(gsub(pattern = "\n|\t",replacement = ",",my_npmi),","))
      my_npmi_spl <- my_npmi_vc[3:length(my_npmi_vc)]
      my_npmi_dt <- data.frame(geneFunc=my_npmi_spl[is.na(as.numeric(my_npmi_spl))],
                               npmi=my_npmi_spl[!is.na(as.numeric(my_npmi_spl))])
      my_npmi_dt$npmi <- as.numeric(my_npmi_dt$npmi)
      my_npmi_dt %>%
        filter(npmi > 0) %>%
        e_charts(geneFunc) %>%
        e_pie(npmi,roseType = "radius",legend = F) %>%
        e_legend(center=0, bottom = 0) %>%
        e_title(text = "Cancer Hallmarks Analytics",my_npmi_vc[1],
                textStyle=list(fontSize = 30),
                left="center",
                subtextStyle=list(fontSize=25)) %>%
        e_tooltip(trigger = "item")
    }
    
    myFeature <- eventReactive(input$myGene_search,{
      return(input$myGene)
    })
    
    output$ft_plot <- renderPlot({
      input$myGene_search
      feature <- isolate(toupper(myFeature()))
      return(FeaturePlot(scRNA(),features = feature))
    },height = 600,width = 770)
    
    output$boder_Vln <- renderPlot({
      input$myGene_search
      feature <- isolate(toupper(myFeature()))
      if (feature != "") {
        return(StackedVlnPlot(scRNA(), feature, pt.size=0,cols=my36colors))
      }
    },height = 200,width = 750) s
    
    output$geneInfo <- shiny::renderUI({
      input$myGene_search
      feature <- isolate(toupper(myFeature()))
      if (feature != "") {
        HTML(NCBI_geninfo(feature))
      }
    })
    
    output$plot_npmi <- echarts4r::renderEcharts4r({
      input$myGene_search
      feature <- isolate(toupper(myFeature()))
      if (feature != "") {
        cancer_npmi(feature)
      }
    })
    
    ################################### PAGE - 4 #######################################################
    geneCardsLink <- function(val,name){
      sprintf('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target="_blank" class="btn btn-primary btn-lg btn-block">%s</a>',val,name)
    }
    output$DEG <- DT::renderDataTable({
      if (length(scRNA()) != 0) {
        myEC <- scRNA()
        variableMarker <- FindAllMarkers(object = myEC,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
        variableMarker$geneCards <- geneCardsLink(rownames(variableMarker),rownames(variableMarker))
        variableMarker %>% 
          dplyr::select(-gene) %>% 
          DT::datatable(
            escape = F,
            colnames = c("gene" = 1),
            filter = 'top', 
            options = list(
              pageLength = 10,
              autoWidth = F
            )
            )
      }
    })
    
    
    
    
}

  

