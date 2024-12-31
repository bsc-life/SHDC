library(shiny)
library(igraph)
library(magrittr)
library(visNetwork)
library(data.table)
library(DT)
library(shinydashboard)
library("shinydashboardPlus")
library(shinythemes)


set.seed(5)


## Options to be selected ##
thenameicd<-c("Name","ICD")
thesex<-c("Women","Men")
similaritymetrics<-c("Pearson","Spearman","Cosine Similarity","Euclidean Distance","Manhattan Distance","Canberra Distance")
geneselections<-c("All genes","Union (sDEGs)","Intersection (sDEGs)")
thesign<-c("Both","Positive","Negative")
onlycomorbid<-c("-","Yes","No")
pathsel<-c("Autophagy","Cell Cycle","Cell-Cell communication","Cellular responses to stimuli","Chromatin organization","Circadian Clock","Developmental Biology","Disease",
           "DNA Repair","DNA Replication","Drug ADME","Extracellular matrix organization","Gene expression Transcription","Hemostasis","Immune System","Metabolism",
           "Metabolism of proteins","Metabolism of RNA","Mitochondria All","Muscle contraction","Neuronal System","Organelle biogenesis and maintenance","Programmed Cell Death",
           "Protein localization","Reproduction","Sensory Perception","Signal Transduction","Transport of small molecules","Vesicle-mediated transport")
comorbidities<-fread("pathways/comorbidities.txt",stringsAsFactors = F,sep="\t",header = F)$V1
pathwayparents<-fread("nodes/Pathway_parents.txt",stringsAsFactors = F,sep="\t")
theparents<-c("all",sort(unique(pathwayparents$parents)))
disdrugterms<-fread("nodes/Disease_drug_identifier.txt",stringsAsFactors = F,sep="\t")
thedrugs<-c("all",disdrugterms$term[which(disdrugterms$category=="drug")])
thedis4drug<-c("all",disdrugterms$term[which(disdrugterms$category=="disease")])

# input<-list("nameoricdm"="Name","sexm"="Women","metricm"="Euclidean Distance","geneselm"="All genes","signom"="Both","comorbidm"="Yes","strengthm"=0,"pathway"="Immune System")

ui <- navbarPage("SHDC",theme=shinytheme("flatly"),
                 tabPanel("Complete network",
                          sidebarLayout(
                            sidebarPanel(
                              width = 2,
                              selectInput("nameoricd","Show disease names or ICD10 codes:",thenameicd,selected = thenameicd[1]),
                              selectInput("sex","Select the sex of interest:",thesex,selected = thesex[1]),
                              selectInput("metric","Select the similarity metric of interest:",similaritymetrics,selected = similaritymetrics[4]),
                              selectInput("genesel","Select the list of genes:",geneselections,selected = geneselections[1]),
                              selectInput("signo","Select the sign of the association:",thesign,selected = thesign[1]),
                              sliderInput("strength","Select the similarity threshold:", min=0, max=1, value=0),
                              selectInput("comorbid","Select comorbid transcriptomic similarities:",onlycomorbid,selected = onlycomorbid[2])
                            ),
                            mainPanel(
                              width=10,
                              h2("Disease Transcriptomic Similarity Network", align = "center"),
                              fluidRow(
                                htmlOutput("txt"),
                                box(
                                  textOutput('textnode1'),
                                  textOutput('textneg1'),
                                  tags$head(tags$style("#textneg1{color: #B64539}")),
                                  textOutput('textpos1'),
                                  tags$head(tags$style("#textpos1{color: #4A74A9}")),width=12,align="center"),
                                box(visNetworkOutput("mononet10_plot", height = "700px"),width = 10),
                                box(imageOutput("icdcat",width="100%"),width = 2,align = "center"),
                                DTOutput("mononet10_table",width='95%'),
                                align="center"
                              ) 
                            )
                          )
                 ),
                 tabPanel("Pathway-based networks",
                          sidebarLayout(
                            sidebarPanel(
                              width = 2,
                              selectInput("nameoricdm","Show disease names or ICD10 codes:",thenameicd,selected = thenameicd[1]),
                              selectInput("sexm","Select the sex of interest:",thesex,selected = thesex[1]),
                              selectInput("metricm","Select the similarity metric of interest:",similaritymetrics,selected = similaritymetrics[3]),
                              selectInput("geneselm","Select the list of genes:",geneselections,selected = geneselections[1]),
                              selectInput("pathway","Select the pathway of interest:",pathsel,selected = pathsel[15]),
                              selectInput("signom","Select the sign of the association:",thesign,selected = thesign[1]),
                              sliderInput("strengthm","Select the similarity threshold:", min=0, max=1, value=0),
                              selectInput("comorbidm","Select comorbid transcriptomic similarities:",onlycomorbid,selected = onlycomorbid[2])
                            ),
                            mainPanel(
                              width=10,
                              h2("Disease Transcriptomic Similarity Network by pathway category", align = "center"),
                              fluidRow(
                                htmlOutput("txt"),
                                box(
                                  textOutput('textnode2'),
                                  textOutput('textneg2'),
                                  tags$head(tags$style("#textneg2{color: #B64539}")),
                                  textOutput('textpos2'),
                                  tags$head(tags$style("#textpos2{color: #4A74A9}")),width=12,align="center"),
                                box(visNetworkOutput("multinet10_plot", height = "700px"),width = 10),
                                box(imageOutput("icdcat2",width="100%"),width = 2,align = "center"),
                                DTOutput("multinet10_table",width='95%'),
                                align="center"
                              )
                            )
                          )
                 ),
                 tabPanel("Pathways in comorbidities observed in both sexes",
                          sidebarLayout(
                            sidebarPanel(
                              width = 2,
                              selectInput("pnameoricd","Show disease names or ICD10 codes:",thenameicd,selected = thenameicd[1]),
                              selectInput("pcomor","Select the comorbidity of interest:",comorbidities,selected = comorbidities[11]),
                              selectInput("pparents","Select the pathway category of interest:",theparents,selected = "all"),
                              selectInput("psex","Select the sex of interest:",choices = c("women","men","common","all"),selected = "all"),
                              selectInput("psigno","Select the sign of the association:",choices = c("up","down","all"),selected = "all")
                            ),
                            mainPanel(
                              width=10,
                              h2("Pathways associated to comorbid diseases in women and men", align = "center"),
                              fluidRow(
                                htmlOutput("txt"),
                                box(
                                  textOutput('textnode3'),
                                  textOutput('textwom3'),
                                  tags$head(tags$style("#textwom3{color: #B64539}")),
                                  textOutput('textmen3'),
                                  tags$head(tags$style("#textmen3{color: #4A74A9}")),
                                  textOutput('textcom3'),
                                  tags$head(tags$style("#textcom3{color: #19787F}")),width=12,align="center"),
                                box(visNetworkOutput("pathnet_plot", height = "700px"),width = 10),
                                box(imageOutput("reactomcat",width="100%"),width = 2,align = "center"),
                                DTOutput("pathnet_table",width='95%'),
                                align="center"
                              )
                            )
                          )
                 ),
                 tabPanel("Drug-disease associations",
                          sidebarLayout(
                            sidebarPanel(
                              width = 2,
                              selectInput("dnameoricd","Show disease names or ICD10 codes:",thenameicd,selected = thenameicd[1]),
                              selectInput("dsex","Select the sex of interest:",choices = c("women","men","all"),selected = "all"),
                              selectInput("ddis","Select the disease of interest:",thedis4drug,selected = thedis4drug[1]),
                              selectInput("ddrug","Select the drug of interest:",thedrugs,selected = thedrugs[1])
                            ),
                            mainPanel(
                              width=10,
                              h2("Diseases connected by drugs in women and men", align = "center"),
                              fluidRow(
                                htmlOutput("txt"),
                                box(
                                  textOutput('textnode4'),
                                  textOutput('textcom4'),
                                  textOutput('textneg4'),
                                  tags$head(tags$style("#textneg4{color: #B64539}")),
                                  textOutput('textpos4'),
                                  tags$head(tags$style("#textpos4{color: #4A74A9}")),width=12,align="center"),
                                
                                box(visNetworkOutput("drugnet_plot", height = "700px"),width = 10),
                                box(imageOutput("icdcat3",width="100%"),width = 2,align = "center"),
                                DTOutput("drugnet_table",width='95%'),
                                align="center"
                              )
                            )
                          )
                 ),
                 tabPanel("Documentation",
                          sidebarLayout(
                            sidebarPanel(
                              width = 0
                            ),
                            mainPanel(
                              width=12,
                              h1("Methods:", align = "center"),
                              br(),
                              br(),
                              h3("Network construction:",align="left"),
                              h5("Transcriptional similarities were calculated on the complete lists of the annotated genes, the union of the annotated significantly differentially expressed genes (sDEG), and their intersection based 
                              on differential expression values (logFC). Six similarity metrics were calculated: Pearson’s and Spearman’s coefficients, cosine similarity, and Euclidean, Canberra, and Manhattan distances. For the cosine 
                              similarity and the Euclidean, Canberra, and Manhattan distances empirical p-values were calculated through 10,000 permutations. In brief, for each gene selection (complete list and union and intersection of 
                              sDEGs) a suffling of the logFC values was performed, and the similarities between diseases were calculated. We corrected for multiple testing by Bonferroni approach, and considered as significant those 
                              similarities with an FDR<=0.05. In the case of Euclidean, Canberra, and Manhattan distances, the mean of the random distances was compared with the actual distances, obtaining positive (negative) values 
                              indicating a greater (lesser) similarity than expected by chance. The similarity values – obtained from the comparison between real and random distances in the case of Euclidean, Canberra and Manhattan 
                              distances, and from the coefficients in the case of Pearson and Spearman correlations and cosine similarity – were binarized, converting those coefficients greater than 0 to +1 and those less than 0 to -1.
                              Going one step further, we generated disease networks with the metrics mentioned above for each comparison (sex) using the genes of each Reactome category (29 in total) and the genes associated with mitochondrial processes (extracted 
                              from MitoCarta [PMID:33174596]) separately. In this way, we generated a multilayer network of diseases, where each layer represents the similarity between diseases based on each of the Reactome categories and mitochondrial genes. Thus, 
                              in each generated network, nodes represent diseases, and edges represent similarities between them based on each metric and gene listing.",align="left"),
                              h3("Overlap with epidemiology:",align="left"),
                              h5("To identify the comorbidity relationships recovered by the disease transcriptomic similarity networks (DTSN) generated by comparing similarities between differential expression profiles and their ability 
                                 to explain comorbidity relationships, we made use of previously published epidemiological network [PMID:30737381] (Supplementary Note 1). The overlap between networks was performed on the shared set of 
                                 diseases (present in both the DTSNs and epidemiological networks). Specifically, the overlap of positive and negative transcriptomic similarities with the epidemiological networks was analyzed separately. 
                                 Overlaps were measured by sex (women vs. women, men vs. men, and adjusted vs. adjusted). The significance of the overlap was assessed by Fisher’s tests and randomizations (generating 10,000 random networks 
                                 shuffling the edges of the DTSNs while maintaining the degree distribution).",align="left"),
                              h3("Disease-drug associations:",align="left"),
                              h5("To study the potential sex-specific role of drugs in comorbidities, we retrieved drug targets from the DrugBank [PMID:29126136]. Since the number of targets per drug is relatively small for enrichment analyses, 
                                 we used the protein-protein interaction network extracted from IID [PMID:34755877] - selecting only those protein-protein interactions in humans that have been experimentally verified – to expand the number of 
                                 targets associated to a drug by mapping the targets on the network and selecting the first neighbors of the targets for each drug. We then conducted a GSEA enrichment analysis [PMID:16199517] to associate drugs 
                                 targeting the products of up- or down-regulated genes with the corresponding disease, separately by sex. Disease-drug associations were extracted from the SIDER database [PMID:26481350]. Disease names were 
                                 transformed into ICD10 codes using the Unified Medical Language System [PMID:14681409] and DrugBank IDs were mapped into drug names.",align="left"),
                              br(),
                              br()
                            )
                          )
                 )
)

server <- function(input, output,session) {
  #### Section one - Monolayer ICD10 ####
  ##   @@ @@ @@ @@ @@ @@ @@ @@ @@   ##
  ## Network ##
  output$icdcat <- renderImage({
    filename <- normalizePath(file.path('./www','ICDcategories.png'))
    list(src = filename,align="center",width="100%",height="100%")
  }, deleteFile = FALSE)
  output$mononet10_plot <- renderVisNetwork({
    ## Select the color of the nodes to be displayed in the network
    nodeinf<-fread("nodes/Disease_colors_icd10.txt",stringsAsFactors = F,sep="\t",header=T)
    ## Load the table ##
    if(input$nameoricd=="Name"){usenames<-"DiseaseName" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10_our_names}
    if(input$nameoricd=="ICD"){usenames<-"ICD10code" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10}
    df<-fread(paste("monolayers/",input$sex,"_",usenames,"_",gsub(" ","",input$metric),"_",substr(input$genesel,1,3),".txt",sep=""),stringsAsFactors = F,sep="\t",header = T)
    ## Select the sign to be displaid ##
    if(input$signo=="Both"){df2<-df}
    if(input$signo=="Positive"){df2<-df[which(df$weight>0)]}
    if(input$signo=="Negative"){df2<-df[which(df$weight<0)]}
    ## Select comorbid conditions ##
    if(input$comorbid=="-"){df3<-df2}
    if(input$comorbid!="-"){df3<-df2[which(df2$comorbid==tolower(input$comorbid))]}
    ## Select the minimum similarity
    df3 <- df3[which(abs(df3$weight) >= input$strength), ]
    ## convert table into a graph
    graph <- graph_from_data_frame(df3, directed=FALSE)
    # is_weighted(graph)
    E(graph)$weight <- df3$weight
    E(graph)$color<-df3$colors
    V(graph)$color<-as.character(nodecolor[names(V(graph))])
    # Visualize the communities
    nodes <- data.frame(id = V(graph)$name, title = V(graph)$name, color = V(graph)$color)
    # nodes <- data.frame(id = V(graph)$name, title = V(graph)$name, color = V(graph)$color,size=V(graph)$size)
    nodes <- nodes[order(nodes$id, decreasing = F),]
    edges <- get.data.frame(graph, what="edges")[1:2]
    edges$color <- df3$colors
    edges$value <- abs(df3$weight)
 
    visNetwork(nodes, edges) %>%
      visExport() %>%
      visOptions(highlightNearest = list(enabled=TRUE, degree=1,algorithm="hierarchical",labelOnly=FALSE), 
                 nodesIdSelection = list(enabled=TRUE,style="width: 300px; height: 26px",main="Select your disease of interest")) %>%
      visIgraphLayout() %>%
      visInteraction(multiselect = T) %>%
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
})
  observeEvent(input$current_node_id, {
    visNetworkProxy("mononet10_plot") %>%
      visGetSelectedNodes()
  })
  ## Text ##
  output$textnode1 <- renderText({"Each node represents a disease, colored based on the category they belong to"})
  output$textneg1 <- renderText({"Red edges represent significant negative similarity values"})
  output$textpos1 <- renderText({"Blue edges represent significant positive similarity values"})
  # Table ##
  output$mononet10_table <- renderDT({
    if(input$nameoricd=="Name"){usenames<-"DiseaseName"}
    if(input$nameoricd=="ICD"){usenames<-"ICD10code"}
    df<-fread(paste("monolayers/",input$sex,"_",usenames,"_",gsub(" ","",input$metric),"_",substr(input$genesel,1,3),".txt",sep=""),stringsAsFactors = F,sep="\t",header = T)
    ## Select the sign to be displaid ##
    if(input$signo=="Both"){df2<-df}
    if(input$signo=="Positive"){df2<-df[which(df$weight>0)]}
    if(input$signo=="Negative"){df2<-df[which(df$weight<0)]}
    ## Select comorbid conditions ##
    if(input$comorbid=="-"){df3<-df2}
    if(input$comorbid!="-"){df3<-df2[which(df2$comorbid==tolower(input$comorbid))]}
    ## Select the minimum similarity
    df3 <- df3[which(abs(df3$weight) >= input$strength), ]
    if(is.null(input$mononet10_plot_selected) | input$mononet10_plot_selected == ''){
      colnames(df3) <- c("Disease1","Disease2","weight","comorbid","colors")
      df3
    }else{
      filtered <- df3[(df3$Disease1 == input$mononet10_plot_selected | df3$Disease2 == input$mononet10_plot_selected),]
      colnames(df3) <- c("Disease1","Disease2","weight","comorbid","colors")
      filtered
    }
  })
  #### Section two - Multilayer ICD10 ####
  ##   @@ @@ @@ @@ @@ @@ @@ @@ @@   ##
  ## Network ##
  output$icdcat2 <- renderImage({
    filename <- normalizePath(file.path('./www','ICDcategories.png'))
    list(src = filename,align="center",width="100%",height="100%")
  }, deleteFile = FALSE)
  output$multinet10_plot <- renderVisNetwork({
    ## Select the color of the nodes to be displayed in the network
    nodeinf<-fread("nodes/Disease_colors_icd10.txt",stringsAsFactors = F,sep="\t",header=T)
    ## Load the table ##
    if(input$nameoricdm=="Name"){usenames<-"DiseaseName" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10_our_names}
    if(input$nameoricdm=="ICD"){usenames<-"ICD10code" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10}
    df<-fread(paste("multilayers/",input$sexm,"_",usenames,"_",gsub(" ","",input$metricm),"_",substr(input$geneselm,1,3),".txt",sep=""),stringsAsFactors = F,sep="\t",header = T)
    ## Select the sign to be displaid ##
    if(input$signom=="Both"){df2<-df}
    if(input$signom=="Positive"){df2<-df[which(df$weight>0)]}
    if(input$signom=="Negative"){df2<-df[which(df$weight<0)]}
    ## Select comorbid conditions ##
    if(input$comorbidm=="-"){df3<-df2}
    if(input$comorbidm!="-"){df3<-df2[which(df2$comorbid==tolower(input$comorbidm))]}
    ## Select the pathway of interest ##
    df3<-df3[which(df3$pathway==input$pathway)]
    ## Select the minimum similarity
    df3 <- df3[which(abs(df3$weight) >= input$strengthm), ]
    ## convert table into a graph
    graph <- graph_from_data_frame(df3, directed=FALSE)
    # is_weighted(graph)
    E(graph)$weight <- df3$weight
    E(graph)$color<-df3$colors
    V(graph)$color<-as.character(nodecolor[names(V(graph))])
    # Visualize the communities
    nodes <- data.frame(id = V(graph)$name, title = V(graph)$name, color = V(graph)$color)
    # nodes <- data.frame(id = V(graph)$name, title = V(graph)$name, color = V(graph)$color,size=V(graph)$size)
    nodes <- nodes[order(nodes$id, decreasing = F),]
    edges <- get.data.frame(graph, what="edges")[1:2]
    edges$color <- df3$colors
    edges$value <- abs(df3$weight)
    
    visNetwork(nodes, edges) %>%
      visExport() %>%
      visOptions(highlightNearest = list(enabled=TRUE, degree=1,algorithm="hierarchical",labelOnly=FALSE), 
                 nodesIdSelection = list(enabled=TRUE,style="width: 300px; height: 26px",main="Select your disease of interest")) %>%
      visIgraphLayout() %>%
      visInteraction(multiselect = T) %>%
      visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
  })
  observeEvent(input$current_node_id, {
    visNetworkProxy("multinet10_plot") %>%
      visGetSelectedNodes()
  })
  ## Text ##
  output$textnode2 <- renderText({"Each node represents a disease, colored based on the category they belong to"})
  output$textneg2 <- renderText({"Red edges represent significant negative similarity values, calculated focused on the genes associated to the selected Reactome pathway category"})
  output$textpos2 <- renderText({"Blue edges represent significant positive similarity values, calculated focused on the genes associated to the selected Reactome pathway category"})
  # Table ##
  output$multinet10_table <- renderDT({
    if(input$nameoricdm=="Name"){usenames<-"DiseaseName"}
    if(input$nameoricdm=="ICD"){usenames<-"ICD10code"}
    df<-fread(paste("multilayers/",input$sexm,"_",usenames,"_",gsub(" ","",input$metricm),"_",substr(input$geneselm,1,3),".txt",sep=""),stringsAsFactors = F,sep="\t",header = T)
    ## Select the sign to be displaid ##
    if(input$signom=="Both"){df2<-df}
    if(input$signom=="Positive"){df2<-df[which(df$weight>0)]}
    if(input$signom=="Negative"){df2<-df[which(df$weight<0)]}
    ## Select comorbid conditions ##
    if(input$comorbidm=="-"){df3<-df2}
    if(input$comorbidm!="-"){df3<-df2[which(df2$comorbid==tolower(input$comorbidm))]}
    ## Select the pathway of interest ##
    df3<-df3[which(df3$pathway==input$pathway)]
    ## Select the minimum similarity
    df3 <- df3[which(abs(df3$weight) >= input$strengthm), ]
    if(is.null(input$multinet10_plot_selected) | input$multinet10_plot_selected == ''){
      colnames(df3) <- c("Disease1","Disease2","weight","comorbid","colors","pathway")
      df3
    }else{
      filtered <- df3[(df3$Disease1 == input$multinet10_plot_selected | df3$Disease2 == input$multinet10_plot_selected),]
      colnames(df3) <- c("Disease1","Disease2","weight","comorbid","colors","pathway")
      filtered
    }
  })
  
  #### Section three - Pathways involved in comorbidities ####
  ## @ @ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ @ @ ##
  # input<-list("dnameoricd"="Name","dsex"="all","ddis"="all","ddrug"="all")
  output$reactomcat <- renderImage({
    filename <- normalizePath(file.path('./www','ReactomeParents.png'))
    list(src = filename,align="center",width="100%",height="100%")
  }, deleteFile = FALSE)
  ## Network ##
  output$pathnet_plot <- renderVisNetwork({
    ## Which is the comorbidity to be displayed? ##
    identifier<-which(comorbidities==input$pcomor)
    ## Select the color of the nodes to be displayed in the network
    nodeinf<-fread("nodes/Disease_colors_icd10.txt",stringsAsFactors = F,sep="\t",header=T)
    ## Load the table ##
    if(input$pnameoricd=="Name"){usenames<-"DiseaseName" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10_our_names}
    if(input$pnameoricd=="ICD"){usenames<-"ICD10code" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10}
    pathcols<-fread("nodes/Pathway_colors.txt",stringsAsFactors = F,sep="\t")
    colpath<-pathcols$colors ; names(colpath)<-pathcols$pathway
    nodecolor<-c(nodecolor,colpath)
    df<-fread(paste("pathways/",usenames,"_network.txt",sep=""),stringsAsFactors = F,sep="\t",header = T)
    ## Select the comorbidity of interest ##
    df<-df[which(df$id==identifier)]
    ## Select the sex of interest ##
    if(input$psex=="all"){df2<-df}
    if(input$psex=="common"){df2<-df[which(df$sex=="common")]}
    if(input$psex=="men"){df2<-df[which(df$sex=="men")]}
    if(input$psex=="women"){df2<-df[which(df$sex=="women")]}
    ## Select the direction of interest ##
    if(input$psigno=="all"){df3<-df2}
    if(input$psigno=="up"){df3<-df2[which(df2$alteration=="up")]}
    if(input$psigno=="down"){df3<-df2[which(df2$alteration=="down")]}
    pathstoselect<-pathwayparents$pathway[which(pathwayparents$parents==input$pparents)]
    if(input$pparents!="all"){
      setkey(df3,"pathway")
      df3<-df3[intersect(df3$pathway,pathstoselect)]
    }
    if(dim(df3)[1]>0){
      ## convert table into a graph
      graph <- graph_from_data_frame(df3, directed=FALSE)
      E(graph)$color<-df3$colors
      colornode<-as.character(nodecolor[names(V(graph))])
      # if(length(which(is.na(colornode)))>0){colornode[which(is.na(colornode))]<-"black"}
      V(graph)$color<-colornode
      V(graph)$shape<-c(rep("circle",2),rep("square",length(as_ids(V(graph)))-2))
      # Visualize the communities
      nodes <- data.frame(id = V(graph)$name, label = V(graph)$name, color = V(graph)$color,shape=V(graph)$shape)
      nodes <- nodes[order(nodes$id, decreasing = F),]
      edges <- get.data.frame(graph, what="edges")[1:2]
      edges$color <- df3$colors
      # edges$lty <- eshapes
      
      visNetwork(nodes, edges) %>%
        visExport() %>%
        visOptions(highlightNearest = list(enabled=TRUE, degree=1,algorithm="hierarchical",labelOnly=FALSE), 
                   nodesIdSelection = list(enabled=TRUE,style="width: 300px; height: 26px",main="Select your node of interest")) %>%
        visIgraphLayout() %>%
        visInteraction(multiselect = T) %>%
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}") 
    }
  })
  observeEvent(input$current_node_id, {
    visNetworkProxy("pathnet_plot") %>%
      visGetSelectedNodes()
  })
  ## Text ##
  output$textnode3 <- renderText({"Circles denote diseases and squares denote pathways"})
  output$textwom3 <- renderText({"Red edges connect pathways to diseases in which they are significantly altered in the same direction only in women"})
  output$textmen3 <- renderText({"Blue edges connect pathways to diseases in which they are significantly altered in the same direction only in men"})
  output$textcom3 <- renderText({"Green edges connect pathways to diseases in which they are significantly altered in the same direction in both women and men"})
  # Table ##
  output$pathnet_table <- renderDT({
    ## Which is the comorbidity to be displayed? ##
    identifier<-which(comorbidities==input$pcomor)
    ## Load the table ##
    if(input$pnameoricd=="Name"){usenames<-"DiseaseName"}
    if(input$pnameoricd=="ICD"){usenames<-"ICD10code"}
    df<-fread(paste("pathways/",usenames,"_network.txt",sep=""),stringsAsFactors = F,sep="\t",header = T)
    ## Select the comorbidity of interest ##
    df<-df[which(df$id==identifier)]
    ## Select the sex of interest ##
    if(input$psex=="all"){df2<-df}
    if(input$psex=="common"){df2<-df[which(df$sex=="common")]}
    if(input$psex=="men"){df2<-df[which(df$sex=="men")]}
    if(input$psex=="women"){df2<-df[which(df$sex=="women")]}
    ## Select the direction of interest ##
    if(input$psigno=="all"){df3<-df2}
    if(input$psigno=="up"){df3<-df2[which(df2$alteration=="up")]}
    if(input$psigno=="down"){df3<-df2[which(df2$alteration=="down")]}
    pathstoselect<-pathwayparents$pathway[which(pathwayparents$parents==input$pparents)]
    if(input$pparents!="all"){
      setkey(df3,"pathway")
      df3<-df3[intersect(df3$pathway,pathstoselect)]
    }
    
    if(is.null(input$pathnet_plot_selected) | input$pathnet_plot_selected == ''){
      colnames(df3) <- c("Disease","pathway","sex","alteration","colors","id")
      df3
    }else{
      filtered <- df3[(df3$Disease == input$pathnet_plot_selected | df3$pathway == input$pathnet_plot_selected),]
      colnames(df3) <- c("Disease","pathway","sex","alteration","colors","id")
      filtered
    }
  })
  #### Section four - Disease - drug associations ####
  ## @@ @@ @@ @@ @@ @@ @ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
  # input<-list("dnameoricd"="Name","dsex"="all","ddis"="E10 (Type 1 diabetes mellitus)","ddrug"="salbutamol")
  output$icdcat3 <- renderImage({
    filename <- normalizePath(file.path('./www','ICDcategories.png'))
    list(src = filename,align="center",width="100%",height="100%")
  }, deleteFile = FALSE)
  ## Network ##
  output$drugnet_plot <- renderVisNetwork({
    ## Select the color of the nodes to be displayed in the network
    nodeinf<-fread("nodes/Disease_colors_icd10.txt",stringsAsFactors = F,sep="\t",header=T)
    ## load the disease-drug associations ##
    disid<-fread("nodes/Disease_id_for_drugs.txt",stringsAsFactors = F,sep="\t")
    ## Load the table ##
    if(input$dnameoricd=="Name"){
      usenames<-"DiseaseName" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10_our_names
      theids<-disid$id ; names(theids)<-disid$diseasename
      if(input$ddis!="all"){thedisease<-as.numeric(theids[gsub("\\)","",gsub(".+\\(","",input$ddis))])}
    }
    if(input$dnameoricd=="ICD"){
      usenames<-"ICD10code" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10
      theids<-disid$id ; names(theids)<-disid$icd10
      if(input$ddis!="all"){thedisease<-as.numeric(theids[gsub(" \\(.+","",input$ddis)])}
    }
    ## colors for indicated diseases ##
    indicatedcols<-fread("nodes/Indicated_diseases_colors.txt",stringsAsFactors = F,sep="\t")
    colind<-indicatedcols$colors ; names(colind)<-indicatedcols$disease
    ## colors for drugs ##
    drugcol<-rep("#000000",length(thedrugs)) ; names(drugcol)<-thedrugs
    ## combine the color vectors ##
    nodecolor<-c(nodecolor,colind,drugcol)
    ## Load the table ##
    df<-fread(paste("drugs/",usenames,"_network.txt",sep=""),stringsAsFactors = F,sep="\t",header = T)
    ## Select the sex of interest ##
    if(input$dsex=="all"){df2<-df}
    if(input$dsex!="all"){df2<-df[which(df$sex==input$dsex)]}
    ## Select the disease of interest ##
    if(input$ddis=="all"){df3<-df2}
    if(input$ddis!="all"){df3<-df2[which(df2$ids==thedisease)]}
    ## Select the drug of interest ##
    if(input$ddrug=="all"){df4<-df3}
    if(input$ddrug!="all"){df4<-df3[unique(c(which(df3$term1==input$ddrug),which(df3$term2==input$ddrug)))]}
    
    if(dim(df4)[1]>0){
      ## convert table into a graph
      graph <- graph_from_data_frame(df4, directed=TRUE)
      E(graph)$weight <- df4$nes
      E(graph)$color<-df4$colors
      colornode<-as.character(nodecolor[names(V(graph))])
      V(graph)$color<-colornode
      nodes <- data.frame(id = V(graph)$name, label = V(graph)$name, color = V(graph)$color)
      nodes <- nodes[order(nodes$id, decreasing = F),]
      edges <- get.data.frame(graph, what="edges")[1:2]
      edges$color <- df4$colors
      edges$weight <- df4$nes
      visNetwork(nodes, edges) %>%
        visExport() %>%
        visEdges(arrows ="to") %>%
        visOptions(highlightNearest = list(enabled=TRUE, degree=1,algorithm="hierarchical",labelOnly=FALSE), 
                   nodesIdSelection = list(enabled=TRUE,style="width: 300px; height: 26px",main="Select your node of interest")) %>%
        visIgraphLayout() %>%
        visInteraction(multiselect = T) %>%
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}") 
    }
  })
  observeEvent(input$current_node_id, {
    visNetworkProxy("drugnet_plot") %>%
      visGetSelectedNodes()
  })
  ## Text ##
  output$textnode4 <- renderText({"Black nodes represent drugs, the other nodes represent diseases"})
  output$textcom4 <- renderText({"Black edges connect diseases to the drugs that treat them (according to SIDER)"})
  output$textneg4 <- renderText({"Red edges connect drugs to the diseases in which they are enriched with a negative normalized enrichment score"})
  output$textpos4 <- renderText({"Blue edges connect drugs to the diseases in which they are enriched with a positive normalized enrichment score"})
  # Table ##
  output$drugnet_table <- renderDT({
    ## Select the color of the nodes to be displayed in the network
    nodeinf<-fread("nodes/Disease_colors_icd10.txt",stringsAsFactors = F,sep="\t",header=T)
    ## load the disease-drug associations ##
    disid<-fread("nodes/Disease_id_for_drugs.txt",stringsAsFactors = F,sep="\t")
    ## Load the table ##
    if(input$dnameoricd=="Name"){
      usenames<-"DiseaseName" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10_our_names
      theids<-disid$id ; names(theids)<-disid$diseasename
      if(input$ddis!="all"){thedisease<-as.numeric(theids[gsub("\\)","",gsub(".+\\(","",input$ddis))])}
    }
    if(input$dnameoricd=="ICD"){
      usenames<-"ICD10code" ; nodecolor<-nodeinf$ICD10_colors ; names(nodecolor)<-nodeinf$ICD10
      theids<-disid$id ; names(theids)<-disid$icd10
      if(input$ddis!="all"){thedisease<-as.numeric(theids[gsub(" \\(.+","",input$ddis)])}
    }
    ## colors for indicated diseases ##
    indicatedcols<-fread("nodes/Indicated_diseases_colors.txt",stringsAsFactors = F,sep="\t")
    colind<-indicatedcols$colors ; names(colind)<-indicatedcols$disease
    ## colors for drugs ##
    drugcol<-rep("#000000",length(thedrugs)) ; names(drugcol)<-thedrugs
    ## combine the color vectors ##
    nodecolor<-c(nodecolor,colind,drugcol)
    ## Load the table ##
    df<-fread(paste("drugs/",usenames,"_network.txt",sep=""),stringsAsFactors = F,sep="\t",header = T)
    ## Select the sex of interest ##
    if(input$dsex=="all"){df2<-df}
    if(input$dsex!="all"){df2<-df[which(df$sex==input$dsex)]}
    ## Select the disease of interest ##
    if(input$ddis=="all"){df3<-df2}
    if(input$ddis!="all"){df3<-df2[which(df2$ids==thedisease)]}
    ## Select the drug of interest ##
    if(input$ddrug=="all"){df4<-df3}
    if(input$ddrug!="all"){df4<-df3[unique(c(which(df3$term1==input$ddrug),which(df3$term2==input$ddrug)))]}
    
    if(is.null(input$drugnet_plot_selected) | input$drugnet_plot_selected == ''){
      colnames(df4) <- c("term1","term2","sex","nes","ids","colors")
      df4
    }else{
      filtered <- df4[(df4$term1 == input$drugnet_plot_selected | df4$term2 == input$drugnet_plot_selected),]
      colnames(df4) <- c("term1","term2","sex","nes","ids","colors")
      filtered
    }
  })
}
shinyApp(ui = ui, server = server)



# To share the library
# library(rsconnect)
# rsconnect::deployApp('/Users/jonsanchezvalle/Desktop/BSC/Griffols/ATshiny')























