library(shiny)
 
shiny::navbarPage(p("scRNA_webGadget",style="font-size:150%;font-weight:bold"),inverse = F,collapsible = F,
                  tabPanel(p("Dat_",align="center",style="font-size:150%;font-weight:bold;margin-top:20%",),
                           fluidRow(style="margin-top:1%;",
                                    column(8,
                                           wellPanel(h1("Summary Cell info",style="margin-top: 0%"),style="background-color: #fff;border:10px double #cccccc; border-bottom: none;height: 260px",
                                                     br(),
                                                     shiny::uiOutput("summary")
                                           ),
                                           wellPanel(h3("Meta.data",style="margin-top: 0%;text-align: center"),
                                                     style="margin-top:-3%;background-color: #fff;border-color: #ffffff;border-left:10px double #cccccc;border-right:10px double #cccccc;",
                                                     br(),
                                                     shinycustomloader::withLoader(DT::dataTableOutput("metadata"))
                                           )
                                    ),
                                    column(4,
                                           fluidRow(
                                               wellPanel(style="background-color: #fff;border:5px solid #cccccc; height: 350px;margin-right:5%",
                                                         shiny::plotOutput("plot_1",width = "100%",height = "100%")
                                               )),
                                           fluidRow(
                                               wellPanel(style="background-color: #fff;border:10px solid #cccccc;height: 450px;margin-right:5%",
                                                         shiny::fileInput("my_file",h3("Select your file")),
                                                         shinyWidgets::switchInput(inputId = "myIQR",label = "Use IQR",labelWidth = "50%",value = T),
                                                         fluidRow(
                                                             helpText("If you turn off the IQR switch above, you will not be able to filter your dataset using the IQR method, but you can customize the filtering criteria by filling in the desired amount of gene expression in the input box.",style="margin-right:40%;margin-left:3%"),
                                                             column(3,numericInput("lower",label = "lower",value = 500,step = 10),style=";margin-left:1%"),
                                                             column(1,p("to",align = "left",style="font-size: 150%;margin-top: 100%;font-weight: bold")),
                                                             column(3,numericInput("upper",label = "upper",step = 1,value = "20000"))
                                                         ),
                                                         br(),
                                                         shiny::actionButton(inputId = "submit","Analysis",style=";margin-left:1%;width:100%")
                                               ))
                                    )
                           )
            ),
            tabPanel(p("Pcs_",align="center",style="font-size:150%;font-weight:bold;margin-top:20%"),
                     fluidRow(style="margin-top:1%;",
                              column(6,
                                     wellPanel(h1("3D-PCA",style="text-align:center;margin-top:1%"),style="background-color: #fff;border:5px outset #cccccc;width:825px ;height: 835px;margin-left: 10%",
                                               threejs::scatterplotThreeOutput(outputId = "plot_pca",width = "100%",height = "90%")
                                               
                                     )
                                     
                              ),
                              
                              column(6,style="margin-left:-5%",
                                     fluidRow(
                                         column(6,
                                                wellPanel(style="background-color: #fff;border:5px dotted #cccccc;width:500px ;height: 400px;margin-left: 10%",
                                                          shinycustomloader::withLoader(shiny::plotOutput("UMAP"))
                                                )
                                                
                                         ),
                                         column(6,style="margin-left:-2%",
                                                wellPanel(style="background-color: #fff;border:5px dotted #cccccc;width:500px ;height: 400px;margin-left: 10%",
                                                          shinycustomloader::withLoader(shiny::plotOutput("tSNE"))
                                                )
                                                
                                         )
                                     ),
                                     fluidRow(
                                         wellPanel(style="border: none;width:95%;margin-left:6%",
                                                   fluidRow(
                                                       shiny::helpText("Here you can customize the principal component dimensions of each axis of the 3D PCA interactive graph."),
                                                       column(4,shinyWidgets::pickerInput("PCA_X",label = "X",choices = paste0("PC_",seq(1,40)),selected = "PC_1" )),
                                                       column(4,shinyWidgets::pickerInput("PCA_Z",label = "Z",choices = paste0("PC_",seq(1,40)),selected = "PC_3" )),
                                                       column(4,shinyWidgets::pickerInput("PCA_Y",label = "Y",choices = paste0("PC_",seq(1,40)),selected = "PC_2" ))
                                                   ),
                                                   fluidRow(
                                                       column(6,
                                                              shinyWidgets::sliderTextInput(
                                                                  inputId = "res",
                                                                  label = h3("Choose a resolution:"),
                                                                  choices = seq(0,2,0.1),
                                                                  selected = 0.5,
                                                                  grid = TRUE,
                                                                  width = '100%',
                                                                  
                                                              )
                                                       ),
                                                       column(6,
                                                              shiny::actionButton("refresh",h3("Refresh",style="text-align:center"),width = "50%",style="margin-left:25%;margin-top:12.5%")
                                                       )
                                                   ),
                                                   fluidRow(
                                                       shinyWidgets::sliderTextInput(
                                                           inputId = "dims",
                                                           label = h3("Choose a dims range to make cluster:"),
                                                           choices = seq(1,40,1),
                                                           selected = c(1,20),
                                                           grid = T,
                                                           
                                                       )
                                                   )
                                                   
                                         )
                                     )
                                     
                              )
                     )
            ),
            tabPanel(p("Fea_",align="center",style="font-size:150%;font-weight:bold;margin-top:20%"),
                     column(6,style="margin-top:1%;",
                            
                            wellPanel(p("Feature Plot",style="text-align:center;margin-top:1%;font-size:200%;font-weight:bolder;color:#fffdfe"),
                                      style="background-color: #fff;border:3px solid #0073b6;width:800px ;height: 800px;margin-left: 9%",
                                      fluidRow(style="margin-top:-9%;",
                                               shiny::plotOutput("boder_Vln"),
                                      ),
                                      fluidRow(style="margin-top:-29%;margin-left:-2%",
                                               shinycustomloader::withLoader(shiny::plotOutput("ft_plot"))
                                      )
                            ),
                            wellPanel(style="background-color: #407ccc;border:1px solid #cccccc;width:800px ;height: 120px;margin-left: 9%",
                                      fluidRow(helpText("Enter the name of the gene you want to query in the input box below, and then click the Go button to draw a map of the gene distribution. At the top, there is a more detailed map of the small request distribution.",
                                                        style="color:#fffdfe")),
                                      fluidRow(
                                          shinyWidgets::searchInput(
                                              inputId = "myGene",
                                              placeholder = "your search gene symbol",
                                              btnSearch = icon("search"), 
                                              btnReset = icon("remove"),
                                              value = "",
                                              width = "80%"
                                          )
                                      ))
                     ),
                     column(5,style="margin-top:1%;margin-left: -4%;",
                            fluidPage(
                                wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:800px ;height: 550px;",
                                          echarts4r::echarts4rOutput(outputId = "plot_npmi",height = "100%")
                                )
                            ),
                            fluidPage(
                                wellPanel(h1("Gene Info",style="text-align:center;font-size:200%;font-weight:bolder;color:#515151;margin-top:-1%"),
                                          style="background-color: #e6e6e6;border:1px solid #cccccc;width:800px ;height: 390px;margin-top:-2%",
                                          shiny::uiOutput(outputId = "geneInfo")
                                )
                            )
                            
                     ),
                     column(1,style="margin-top:1%;margin-left: 1%;",
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;"),
                            wellPanel(style="background-color: #fff;border:1px solid #cccccc;width:1% ;height: 20px;")
                     )
            ),
            tabPanel(p("DEG_",align="center",style="font-size:150%;font-weight:bold;margin-top:20%"),
                     column(6,style="margin-top:1%;",
                            wellPanel(style="width:200%;height:800px",
                                      DT::dataTableOutput("DEG")
                            )
                     )
            )
)



