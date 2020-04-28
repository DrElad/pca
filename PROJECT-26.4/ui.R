library(shiny) 
library(shinyBS)
library(readxl)
library(ggplot2)
library(Cairo)
library(shinyjs)
library(htmlwidgets)

#################UI#######################
ui <- fluidPage(
 
  textInput(inputId = "my_gene", label = "Enter gene name", value = "", width = "300%", placeholder = NULL), #creating a new variable
  textInput(inputId = "PCA_X", label = "Enter PCA in x axis", value = "", width = "300%", placeholder = NULL), #creating a new variable
  textInput(inputId = "PCA_Y", label = "Enter PCA in y axis", value = "", width = "300%", placeholder = NULL), #creating a new variable
  actionButton(inputId = "go", label = "search gene"), #creating a button that shows the relevant gene's row in the "TheDataArrenge"
 ############################################### 


################default######################
  fluidRow(
    column(width = 10,
           plotOutput("plot4", height = 600,
                      # Equivalent to: click = clickOpts(id = "plot_click")
                      click = "plot4_click",
                      brush = brushOpts(
                        id = "plot4_brush"
                      )
           )
    )
  ),
###################################
  fluidRow(
    column(width = 10,
           plotOutput("plot1", height =600 ,
                      # Equivalent to: click = clickOpts(id = "plot_click")
                      click = "plot1_click",
                      brush = brushOpts(
                        id = "plot1_brush"
                      )
           )
    )
  ),
  
  ############################################
  fluidRow(
    column(width = 5,
           h4("Points near click"),
           verbatimTextOutput("click_info")
    ),
    column(width = 5,
           h4("Brushed points"),
           verbatimTextOutput("brush_info")
    )
   ),
  #######################ZOOMIN##################
  fluidRow(
    column(width = 10, 
           h4("Brush and double-click to zoom"),
           plotOutput("plot2", height = 600,
                      dblclick = "plot2_dblclick",
                      brush = brushOpts(
                        id = "plot2_brush"
                        
                      )
           )
    ),
 
    
  ),
  ###################HIST################################
  fluidRow(
    column(width = 10, 
           h4("hist"),
           plotOutput("plot3", height = 500,
                      dblclick = "plot3",
                      brush = brushOpts(
                        id = "plot3"
                        
                      )
           )
    ),
    
    
  )
  ###################SERVER##############################
)
server <- function(input, output) { 
  ########################DATALOAD#########################
  nem_score<-read_excel('./data/score_matrix_with_archetypes.xls',col_name=FALSE)
  TheData<-read_excel('./data/postpartitable_dkfz.xls')
  TheDataArrenge<-TheData[,-1]
  row.names(TheDataArrenge)<-TheData$Row
  ##################creating the default graph############
  Z_def<-TheDataArrenge[c(as.character.default("LRP2")),]
  def_gene<-t(Z_def)
  def_gene<-as.data.frame(def_gene)
  pca_y1_def_gene<-nem_score[,1]
  pca_y3_def_gene<-nem_score[,2]
  plot_def_gene<-data.frame(def_gene,pca_y1_def_gene,pca_y3_def_gene)
  names(plot_def_gene)[1]<-as.character("LRP2")
  names(plot_def_gene)[2]<-"PCA1"
  names(plot_def_gene)[3]<-"PCA2"
  def_gene_exp_level<-plot_def_gene[,1]
  output$plot4<-renderPlot(ggplot(plot_def_gene, aes(x = plot_def_gene[,2], y = plot_def_gene[,3]))+
  geom_point(aes( color=def_gene_exp_level,size =def_gene_exp_level))+ggtitle(as.character("LRP2"))+theme(plot.title = element_text(colour="red", size="14" , face="bold.italic"))+scale_colour_gradient(low = "green", high = "red")+xlab('PCA1')+ylab('PCA2'))
 #################MAIN PLOT#############################
    observeEvent(input$go,{
      n<-input$my_gene
      pca_x<-input$PCA_X
      pca_y<-input$PCA_Y
      Z<-TheDataArrenge[c(as.character(input$my_gene)),]
      pca_x_1<-nem_score[c(as.numeric(input$PCA_X))]
      pca_y_1<-nem_score[c(as.numeric(input$PCA_Y))]
      Z<-TheDataArrenge[c(as.character(n)),]
      pca_x_1<-nem_score[c(as.numeric(pca_x))]
      pca_y_1<-nem_score[c(as.numeric(pca_y))]
      THEGENE<-t(Z)
      THEGENE<-as.data.frame(THEGENE)
      pca_X_1_1<-as.data.frame(pca_x_1)
      pca_y_1_1<-as.data.frame(pca_y_1)
      pca_x_f<-data.frame(pca_X_1_1)
      pca_y_f<-data.frame(pca_y_1_1)
      main_plot<-data.frame(THEGENE,pca_x_f,pca_y_f)
      names(main_plot)[1]<-as.character(n)
      names(main_plot)[2]<- "PCA X"
      names(main_plot)[3]<- "PCA Y"
       gene_level_expresion<-main_plot[,1]
       output$plot1<-renderPlot(ggplot(main_plot,aes(x =main_plot[,2],y = main_plot[,3]))+
        geom_point(aes( color=gene_level_expresion,size=gene_level_expresion))+ggtitle(as.character(n))+
          theme(plot.title = element_text(colour="red", size="14", face="bold.italic"))+scale_colour_gradient(low = "green", high = "red")+
          xlab("PCA X")+ylab("PCA Y"))
      #######################################
      
  #############CLICK_ON_POINT##########
  output$click_info <- renderPrint({
    # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
    # were a base graphics plot, we'd need those.
      nearPoints(main_plot, input$plot1_click,"PCA X","PCA Y")
  })
  ############BRUSH####################
  output$brush_info <- renderPrint({
    brushedPoints(main_plot, input$plot1_brush,"PCA X","PCA Y")
  })
  ###############ZOOMIN#############
  
 ##############ZOOM_IN###################
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$plot2 <- renderPlot({
    ggplot(main_plot, aes(main_plot[,2], main_plot[,3])) +
      geom_point(aes( color=gene_level_expresion,size =gene_level_expresion)) +
      coord_cartesian(xlim = ranges$x, ylim = ranges$y)+scale_colour_gradient(low = "green", high = "red")+xlab(as.numeric(input$PCA_X))+ylab(as.numeric(input$PCA_Y))
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot2_dblclick, {
    brush <- input$plot2_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  #####################HIST############################
  output$plot3 <- renderPlot({
    
      hist(x =as.numeric(Z)  , col = 'blue', border = 'white',main = as.character(n) ,xlab =as.character(n))
  }) 
  })

  } 

shinyApp(ui = ui, server = server,options = list(height = 1080))