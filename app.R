#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(leaflet)
library(sp)
library(shiny)
library(shinyWidgets)
library(geojsonio)

columnas <- c('Species','Gene ID','Complete Gene ID','Gene Name and Gene Symbol','PANTHER family / subname',
  'PANTHER protein class','GO-slim Molecular Function','GO-slim Biological Process','GO-slim Cellular Component',
  'GO database M.F. (complete)','GO database B.P. (complete)','GO database C.C. (complete)')

col_boxes <- append(columnas, c('logFC', 'logCPM', 'PValue', 'FDR'))

upreg_num <- NULL
downreg_num <- NULL

deg_up_3 <- read.delim("upreg/pantherGeneList_1.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[3] <- nrow(deg_up_3)
deg_up_5 <- read.delim("upreg/pantherGeneList_1.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[5] <- nrow(deg_up_5)
deg_up_15 <- read.delim("upreg/pantherGeneList_2.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[15] <- nrow(deg_up_15)
deg_up_17 <- read.delim("upreg/pantherGeneList_3.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[17] <- nrow(deg_up_17)
deg_up_27 <- read.delim("upreg/pantherGeneList_4.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[27] <- nrow(deg_up_27)
deg_up_31 <- read.delim("upreg/pantherGeneList_5.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[31] <- nrow(deg_up_31)
deg_up_35 <- read.delim("upreg/pantherGeneList_6.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[35] <- nrow(deg_up_35)
deg_up_36 <- read.delim("upreg/pantherGeneList_7.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[36] <- nrow(deg_up_35)
deg_up_38 <- read.delim("upreg/pantherGeneList_8.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[38] <- nrow(deg_up_38)
deg_up_41 <- read.delim("upreg/pantherGeneList_9.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
upreg_num[41] <- nrow(deg_up_41)

deg_down_3 <- read.delim("downreg/pantherGeneList_1.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[3]<- nrow(deg_down_3)
deg_down_5 <- read.delim("downreg/pantherGeneList_1.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[5]<- nrow(deg_down_5)
deg_down_15 <- read.delim("downreg/pantherGeneList_2.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[15] <- nrow(deg_down_15)
deg_down_17 <- read.delim("downreg/pantherGeneList_3.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[17] <- nrow(deg_down_17)
deg_down_27 <- read.delim("downreg/pantherGeneList_4.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[27] <- nrow(deg_down_27)
deg_down_31 <- read.delim("downreg/pantherGeneList_5.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[31] <- nrow(deg_down_31)
deg_down_35 <- read.delim("downreg/pantherGeneList_6.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[35] <- nrow(deg_down_35)
deg_down_36 <- read.delim("downreg/pantherGeneList_7.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[36] <- nrow(deg_down_36)
deg_down_38 <- read.delim("downreg/pantherGeneList_8.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[38] <- nrow(deg_down_38)
deg_down_41 <- read.delim("downreg/pantherGeneList_9.txt", sep="\t", header=FALSE, col.names = columnas, check.names = FALSE)
downreg_num[41] <- nrow(deg_down_41)

crop_polygons <- geojson_read("map.geojson",
                      what = "sp")
stages_year <- read.table(file = "stages_year.tsv", sep = "\t")
colnames(stages_year) <-  c("Agosto 2017", "Septiembre 2017", "Octubre 2017", "Noviembre 2017", "Diciembre 2017", 
                        "Enero 2018", "Febrero 2018", "Marzo 2018", "Abril 2018", "Mayo 2018", "Junio 2018", "Julio 2018")
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Crop development analysis using RNA-Seq data"),
  
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      
      sliderTextInput("month",
                      "Month of the year:",
                      choices = c("Agosto 2017", "Septiembre 2017", "Octubre 2017", "Noviembre 2017", "Diciembre 2017",
                                  "Enero 2018", "Febrero 2018", "Marzo 2018", "Abril 2018", "Mayo 2018", "Junio 2018", "Julio 2018"),
                      # choices = c(3,5,15,17,27,31,35,36,38,41),
                      selected = "Agosto 2017"
                      #min = 1,
                      #max = 50,
                      #value = 30)
      ),
      
      hr(),
      
      selectInput("stage2show", label = h4("Select stage to see differentially expressed genes"), 
                  choices = list("EL-3" = 3, "EL-5" = 5, "EL-15" = 15, "EL-17" = 17, "EL-27" = 27, 
                                 "EL-31" = 31, "EL-35" = 35, "EL-36" = 36, "EL-38" = 38, "EL-41" = 41), 
                  selected = 3),
      
      # textOutput("texto"),
      checkboxGroupInput("columns2show", "Columns to show:", 
                         col_boxes, selected = c('Gene ID', 'Gene Name and Gene Symbol','PANTHER family / subname','PANTHER protein class','GO-slim Molecular Function','GO-slim Biological Process','GO-slim Cellular Component', 'logFC', 'PValue'))
    ),
    
    
    leafletOutput(outputId = "map",height = 700, width = 650)
    
    # Show a plot of the generated distribution
    #mainPanel(
    #   plotOutput("distPlot")
    #)
  ),
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Upregulated", DT::dataTableOutput("DEG_up")),
                tabPanel("Downregulated", DT::dataTableOutput("DEG_down"))
    )
    
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # output$distPlot <- renderPlot({
  #    # generate bins based on input$bins from ui.R
  #    x    <- faithful[, 2] 
  #    bins <- seq(min(x), max(x), length.out = input$bins + 1)
  #    
  #    # draw the histogram with the specified number of bins
  #    hist(x, breaks = bins, col = 'darkgray', border = 'white')
  # })
  
  output$map <- renderLeaflet({
    stg = c(3,5,15,17,27,31,35,36,38,41)
    #pal <- colorNumeric("YlOrRd", domain = stg)
    pal <- colorBin(palette = "YlOrRd", domain = stg)

    if(!is.na(stages_year[crop_polygons$zona_n,input$month])){
      labels <- sprintf(
        "<strong>Sector %g: Vitis vinifera</strong><br/>Development Stage: EL-%g<br/>%g Differentially Expressed Genes<br/>Upregulated: %g, Downregulated: %g",
          crop_polygons$zona_n, stages_year[crop_polygons$zona_n,input$month], upreg_num[stages_year[crop_polygons$zona_n,input$month]] + downreg_num[stages_year[crop_polygons$zona_n,input$month]], upreg_num[stages_year[crop_polygons$zona_n,input$month]], downreg_num[stages_year[crop_polygons$zona_n,input$month]]
        ) %>% lapply(htmltools::HTML)
    }
    else{
      labels <- sprintf(
        "Harvested"
        ) %>% lapply(htmltools::HTML)
    }

    m <- leaflet(data = crop_polygons) %>% 
      
      addPolygons(color = "black", fillOpacity = 0.8, highlightOptions = highlightOptions(
                                     color = "white", weight = 2, bringToFront = TRUE),
                  fillColor = pal(stages_year[crop_polygons$zona_n,input$month]),
                  label = labels,
                  labelOptions = labelOptions(
                    style = list("font-weight" = "normal", padding = "3px 8px"),
                    textsize = "15px",
                    direction = "auto")
                  ) %>% 
      
      addProviderTiles(providers$Esri.WorldImagery) %>% 
      addLegend(pal = pal, values = stg, opacity = 0.7, title = "EL Stage", position = "bottomright")
    m
  })
  

  # output$texto <- renderText(
  #   expr <- paste0("seleccionado: ",input$map_shape_click['id'])
  #   )
  
  output$DEG_up <- DT::renderDataTable({
    if(input$stage2show == 3){
      deg_up <- deg_up_3
      #Agregar valores cuantitativos (logFC, logCPM, p-value, FDR)
      up_values <- read.delim("upreg/values/upreg_1.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_3, up_values, by='Gene ID')
    }
    if(input$stage2show == 5){
      deg_up <- deg_up_5
      up_values <- read.delim("upreg/values/upreg_1.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_5, up_values, by='Gene ID')
    }
    if(input$stage2show == 15){
      deg_up <- deg_up_15
      up_values <- read.delim("upreg/values/upreg_2.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_15, up_values, by='Gene ID')
    }
    if(input$stage2show == 17){
      deg_up <- deg_up_17
      up_values <- read.delim("upreg/values/upreg_3.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_17, up_values, by='Gene ID')
    }
    if(input$stage2show == 27){
      deg_up <- deg_up_27
      up_values <- read.delim("upreg/values/upreg_4.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_27, up_values, by='Gene ID')
    }
    if(input$stage2show == 31){
      deg_up <- deg_up_31
      up_values <- read.delim("upreg/values/upreg_5.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_31, up_values, by='Gene ID')
    }
    if(input$stage2show == 35){
      deg_up <- deg_up_35
      up_values <- read.delim("upreg/values/upreg_6.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_35, up_values, by='Gene ID')
    }
    if(input$stage2show == 36){
      deg_up <- deg_up_36
      up_values <- read.delim("upreg/values/upreg_7.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_36, up_values, by='Gene ID')
    }
    if(input$stage2show == 38){
      deg_up <- deg_up_38
      up_values <- read.delim("upreg/values/upreg_8.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_38, up_values, by='Gene ID')
    }
    if(input$stage2show == 41){
      deg_up <- deg_up_41
      up_values <- read.delim("upreg/values/upreg_9.out", sep=" ")
      up_values <- cbind(up_values, 'Gene ID' = rownames(up_values))
      deg_up <- merge(deg_up_41, up_values, by='Gene ID')
    }
    expr <- deg_up
    expr[, input$columns2show, drop = FALSE]
  })

  output$DEG_down <- DT::renderDataTable({
    if(input$stage2show == 3){
      deg_down <- deg_down_3
      down_values <- read.delim("downreg/values/downreg_1.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_3, down_values, by='Gene ID')
    }
    if(input$stage2show == 5){
      deg_down <- deg_down_5
      down_values <- read.delim("downreg/values/downreg_1.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_5, down_values, by='Gene ID')
    }
    if(input$stage2show == 15){
      deg_down <- deg_down_15
      down_values <- read.delim("downreg/values/downreg_2.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_15, down_values, by='Gene ID')
    }
    if(input$stage2show == 17){
      deg_down <- deg_down_17
      down_values <- read.delim("downreg/values/downreg_3.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_17, down_values, by='Gene ID')
    }
    if(input$stage2show == 27){
      deg_down <- deg_down_27
      down_values <- read.delim("downreg/values/downreg_4.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_27, down_values, by='Gene ID')
    }
    if(input$stage2show == 31){
      deg_down <- deg_down_31
      down_values <- read.delim("downreg/values/downreg_5.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_31, down_values, by='Gene ID')
    }
    if(input$stage2show == 35){
      deg_down <- deg_down_35
      down_values <- read.delim("downreg/values/downreg_6.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_35, down_values, by='Gene ID')
    }
    if(input$stage2show == 36){
      deg_down <- deg_down_36
      down_values <- read.delim("downreg/values/downreg_7.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_36, down_values, by='Gene ID')
    }
    if(input$stage2show == 38){
      deg_down <- deg_down_38
      down_values <- read.delim("downreg/values/downreg_8.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_38, down_values, by='Gene ID')
    }
    if(input$stage2show == 41){
      deg_down <- deg_down_41
      down_values <- read.delim("downreg/values/downreg_9.out", sep=" ")
      down_values <- cbind(down_values, 'Gene ID' = rownames(down_values))
      deg_down <- merge(deg_down_41, down_values, by='Gene ID')
    }
    expr <- deg_down
    expr[, input$columns2show, drop = FALSE]
  })
  observeEvent(input$map_shape_click, { # update the location selectInput on map clicks
    p <- input$map_shape_click
    print(p)

  })
}

# Run the application 
shinyApp(ui = ui, server = server)

