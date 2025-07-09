#
# Avispace version 4.0
# 6 June 2025
# Compared to version 3.0, this version is greatly simplified to avoid crashing the public server
# modified for flat polygons, no options to modify eigenvectors

library(shiny)
library(shinythemes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(sf)
library(factoextra)
library(rnaturalearth)
library(rnaturalearthdata)

# working directory (for trial only) -------------------------------------------
#setwd("/Users/lauraberman/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/Documents/Wisconsin/Townsend Lab/Traits and acoustics/Shiny/Avispace4.0")


# read in data -----------------------------------------------------------------
speciesPolygons <- readRDS("speciesPolygons_jittered_detPerDay_20250624.rds")
siteDetections_foliarTraits_BioCube <- readRDS("siteDetections_foliarTraits_BioCube_20250522.rds")
NicheAreas <- readRDS("speciesNicheAreas_20250624.rds")
NicheRangeOverlap <- readRDS("NicheRangeOverlap_dotNames_20250624.rds")
rangeOverlaps <- readRDS("rangeOverlaps_dotNames_20250623.rds")
GeospatialMetadata <- read.csv("GeospatialMetadata.csv")


# convert to count data --------------------------------------------------------
names(siteDetections_foliarTraits_BioCube)
roundedCounts <- siteDetections_foliarTraits_BioCube %>% mutate_at(5:96, round, 0) # double check these are species columns
roundedCounts$geometry <- NULL

# name vars --------------------------------------------------------------------
species <- colnames(roundedCounts)[5:96]
spatVars <- colnames(roundedCounts)[97:182]

# calculate PCA ----------------------------------------------------------------
traitsPCA <-  prcomp(roundedCounts[,spatVars], center=TRUE, scale = TRUE)
# add PCs to rounded dataframe
roundedPCA <- cbind(roundedCounts, traitsPCA$x)

# Function to convert sf POLYGON to data frame for ggplot ----------------------
sf_to_df <- function(sf_obj) {
  sf::st_geometry(sf_obj) %>%
    sf::st_cast("POLYGON") %>%
    lapply(sf::st_coordinates) %>%
    lapply(as.data.frame) %>%
    bind_rows(.id = "group")
}


# UI ###########################################################################

ui <- fluidPage(
  theme = shinytheme("journal"),
  titlePanel(" AviSpace: Avian niche specificity in the Sierra Nevada, USA"),
  h5(" AviSpace is a visualization tool for avian habitat perferences, habitat 
     specificity, and niche overlap, based on data from the Sierra Nevada Mountains 
     in California, USA. Bird detection densities are plotted in PCA niche 
     space. Habitat generalists occur across larger areas of niche space, while 
     habitat specialists occupy smaller areas."), 
  navbarPage("", 
             tabPanel("PCA",
                      sidebarLayout(
                        sidebarPanel( width = 3,
                                      HTML("<label style='color:#25b9be; font-weight:bold;'>Select bird species #1:</label>"),
                          selectInput("Species1", 
                                      NULL,
                                      choices = species, 
                                      selected = "Acorn.Woodpecker"),
                          HTML("<label style='color:#e41a1c; font-weight:bold;'>Select bird species #2:</label>"),
                          selectInput("Species2", 
                                      NULL,
                                      choices = species, 
                                      selected = "Townsend.s.Solitaire"),
                        ),
                        mainPanel(
                          fluidRow(
                            column(width = 7, plotOutput("pcaPlot"),
                                   uiOutput("Species1NicheArea"),
                                   uiOutput("Species2NicheArea"),
                                   uiOutput("NicheOverlap")),
                            column(width = 5, plotOutput("rangeMap"),
                                   uiOutput("Species1RangeArea"),
                                   uiOutput("Species2RangeArea"),
                                   uiOutput("RangeOverlap"))),
                          h6("LEFT: PCA of habitat niche space. Grey points are 
                          individual recording sites. Text indicates general landcover 
         type, more detailed landcover information is available in the 'about' tab. 
         Filled color indicates the minimum area contour where at least 95% of 
         detections of the selected species occur. 
         Bird detections are based on passive acoustic monitoring and identified using BirdNET. 
         Areas are measured in unitless niche space. 
         RIGHT: Breeding range maps of the two selected species, from eBird status and trends.")
                        )
                      )
             ),
             tabPanel("Boxplot",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("Birds",
                                      "Select bird species:",
                                      choices = colnames(roundedCounts)[5:96],
                                      selected = c("Oak.Titmouse", "Townsend.s.Solitaire", "American.Robin", "Golden.crowned.Kinglet"),
                                      multiple = TRUE),
                          selectInput("SpatVars",
                                      "Select spatial variable:",
                                      choices = spatVars,
                                      selected = "Nitrogen")
                        ),
                        mainPanel(
                          plotOutput("boxPlot")
                        ),
                      )),
             tabPanel("Data Download",
                      mainPanel(
                        h5("Niche Areas"),
                        p("The calculated niche area of each species in unitless PCA space."),
                        downloadButton("downloadNicheAreas", "Download Niche Areas"),
                        h5("Niche Range Overlap"),
                        p("The percent niche overlap and range overlap of all pairwise species combinations."),
                        downloadButton("downloadNicheOverlaps", "Download Niche Range Overlap"),
                        h5("Geospatial Metadata"),
                        p("Source details for geospatial layers"),
                        downloadButton("downloadGeospatial", "Download Geospatial Metadata")
                      )
             ),
             tabPanel("About",
                      img(src='methods figure.png', width=800, align = "left"),
                      img(src='landcover_pca_20250522.png', width=800, align = "left"),
                      h5("Study sites"),
                      p("AviSpace is based on data from 577 study sites distributed 
                      throughout the Sierra Nevada mountain range in California, USA."),
                      h5("Bird data"),
                      p("Bird detection densities are based on acoustic detection rates 
                      at autonomous recording units (ARUs) deployed at each of 577 study 
                      sites in 2021-2023. Audio data collection efforts were led by 
                      Connor Wood and Zach Peery."),
                      h5("PCA Spatial Variables"),
                      p("PCA eigenvectors included 92 geospatial variables describing 
                      habitat dimensions including climate, terrain, anthropogenic 
                      disturbance, phenology, leaf physiological traits, and ecosystem function. 
                      Foliar traits were modeled based on hyperspectral imagery 
                      collected in 2018. Trait models were developed primarily by 
                      Ethan Shafron and Ting Zheng. The same PCA ordination was used for all bird species."),
                      h5("Excluded Sites - Fire"),
                      p("Because bird data (2021-2023) and foliar trait data (2018) were collected in different 
                      years, any sites which burned 2018-2023 were excluded."),
                      h5("How to cite AviSpace"),
                      p("If you use AviSpace in your published work, please cite"),
                      p("Berman, LM; Schneider, FD; Pavlick, RP; Peery, MZ; Wood, CM; 
                        Zheng, T; Shafron, E; Ye, Z; Queally, N; Dean, M; Tagliabue, G;
                        Winiarski, J; Kramer, A; Townsend, P. (in prep) 
                        Remote sensing and bioacoustics reveal avian niche partitioning and habitat filtering")
             ),
  )
)


# Server #######################################################################

server <- function(input, output) {
  
  
  # get polygons ---------------------------------------------------------------
  species1_poly <- reactive({speciesPolygons %>% filter(species == input$Species1)}) 
  species2_poly <- reactive({speciesPolygons %>% filter(species == input$Species2)})
  
  # Calculate intersection
  overlap_poly <- reactive({st_intersection(species1_poly(), species2_poly())})
  
  # convert to df
  species1_df <- reactive({sf_to_df(species1_poly())})
  species2_df <- reactive({sf_to_df(species2_poly())})
  overlap_df  <- reactive({sf_to_df(overlap_poly())})
  
  # Plot with overlap ----------------------------------------------------------
  output$pcaPlot <- renderPlot({ 
    
    base_plot <- fviz_pca_ind(traitsPCA,
                 label = "none",
                 alpha.ind = 0.1) +
      geom_polygon(data = species1_df(), aes(x = X, y = Y, group = group), fill = "#25b9be", alpha = 0.6) +
      geom_polygon(data = species2_df(), aes(x = X, y = Y, group = group), fill = "#e41a1c", alpha = 0.6) +
      geom_text(aes(x = -7, y = -8, label = "Oak"), inherit.aes = FALSE) +
      geom_text(aes(x =  6, y = -2, label = "Pine"), inherit.aes = FALSE) +
      geom_text(aes(x =  0, y = 1, label = "Mixed Conifer"), inherit.aes = FALSE) +
      coord_fixed(ratio = 17.2/27.4) +
      theme_minimal() +
      theme(legend.title=element_blank(), legend.position = "bottom") +
      labs(title = "")
    
    # Conditionally add overlap layer if it has data
    if (nrow(overlap_df()) > 0) {
      base_plot <- base_plot +
        geom_polygon(data = overlap_df(), aes(x = X, y = Y, group = group), fill = "#8849a5", alpha = 0.6)
    }
    
    base_plot
    })

  

  # make range maps ------------------------------------------------------------ 
  
  output$rangeMap <- renderPlot({
    
    req(input$Species1, input$Species2)
    
    # find the paths
    map_paths <- list.files("RangeMaps/", 
                            pattern = ".*_range_(2022|2023)\\.gpkg$", full.names = TRUE, recursive = TRUE)
    
    # get the maps
    Species1_map <- st_read(grep(input$Species1, map_paths, value=TRUE))
    Species2_map <- st_read(grep(input$Species2, map_paths, value=TRUE))
    
    # only breeding season please
    Species1_map_breeding <- subset(Species1_map, season %in% c("resident", "breeding"))
    Species2_map_breeding <- subset(Species2_map, season %in% c("resident", "breeding"))
    
    # get intersection
    intersection_map <- st_intersection(Species1_map_breeding, Species2_map_breeding)
    
    # Get North America map
    north_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
    
    # Define a bounding box (xmin, ymin, xmax, ymax)
    bbox <- st_bbox(c(xmin = -140, xmax = -50, ymin = 5, ymax = 70), crs = st_crs(north_america))
    
    # Clip to the bounding box
    north_america_clipped <- st_crop(north_america, bbox)
    Species1_map_clipped <- st_crop(Species1_map_breeding, bbox)
    Species2_map_clipped <- st_crop(Species2_map_breeding, bbox)
    intersection_map_clipped <- st_crop(intersection_map, bbox)
    
    # plot it
    base_map <- ggplot() +
      geom_sf(data = north_america_clipped, fill = "grey40", color = "grey40") +
      geom_sf(data = Species1_map_clipped, fill = "#25b9be", color = "#25b9be", size = 1, alpha=0.6) +
      geom_sf(data = Species2_map_clipped, fill = "#e41a1c", color = "#e41a1c", size = 1, alpha=0.6) +
      theme_void() +  
      theme(legend.position = "none")
    
    # Conditionally add overlap layer if it has data
    if (nrow(intersection_map_clipped) > 0) {
      base_map <- base_map +
        geom_sf(data = intersection_map_clipped, fill = "#8849a5", color = "#8849a5", size = 1, alpha=0.6)
    }
    
    base_map
    
    
  })
    
  
  # get niche areas ------------------------------------------------------------
  
  sp1Area <- reactive({NicheAreas$area[NicheAreas$species == input$Species1]})
  sp2Area <- reactive({NicheAreas$area[NicheAreas$species == input$Species2]})
  sp12NicheOverlap <- reactive({
    val <- NicheRangeOverlap$niche_overlap_percent[NicheRangeOverlap$species1 == input$Species1 & NicheRangeOverlap$species2 == input$Species2]
    if (length(val) == 0) {val <- NicheRangeOverlap$niche_overlap_percent[NicheRangeOverlap$species1 == input$Species2 & NicheRangeOverlap$species2 == input$Species1]}
    val
    })
  
  # output niche areas
  output$Species1NicheArea <- renderUI({HTML(paste0("<span style='color:#25b9be;'>", input$Species1, " Niche Area: ", round(sp1Area(),2), "</span>"))})
  output$Species2NicheArea <- renderUI({HTML(paste0("<span style='color:#e41a1c;'>", input$Species2, " Niche Area: ", round(sp2Area(),2), "</span>"))})
  output$NicheOverlap <- renderUI({HTML(paste0("<span style='color:#8849a5;'>", " Niche Overlap: ", round(sp12NicheOverlap()*100,2), "%", "</span>"))})
  
  
  # get range areas ------------------------------------------------------------
  
  sp1Range <- reactive({rangeOverlaps$area1[rangeOverlaps$species1 == input$Species1][1]})
  sp2Range <- reactive({rangeOverlaps$area1[rangeOverlaps$species1 == input$Species2][1]})
  sp12RangeOverlap <- reactive({
    val <- NicheRangeOverlap$range_overlap_percent[NicheRangeOverlap$species1 == input$Species1 & NicheRangeOverlap$species2 == input$Species2]
    if (length(val) == 0) {val <- NicheRangeOverlap$range_overlap_percent[NicheRangeOverlap$species1 == input$Species2 & NicheRangeOverlap$species2 == input$Species1]}
    val
    })
  
  # output niche areas
  output$Species1RangeArea <- renderUI({HTML(paste0("<span style='color:#25b9be;'>", input$Species1, " Range Area: ", round(sp1Range()/1000000000000,3), " million kM2", "</span>"))})
  output$Species2RangeArea <- renderUI({HTML(paste0("<span style='color:#e41a1c;'>", input$Species2, " Range Area: ", round(sp2Range()/1000000000000,3), " million kM2", "</span>"))})
  output$RangeOverlap <- renderUI({HTML(paste0("<span style='color:#8849a5;'>", " Range Overlap: ", round(sp12RangeOverlap()*100,2), "%", "</span>"))})
  

  
  
  # boxplot panel --------------------------------------------------------------
  # pivot longer for selected species
  selectedSpecies <- reactive({
    pivot_longer(roundedCounts,
                 cols = input$Birds,
                 names_to = "species",
                 values_to = "counts")
  })
  
  # uncount it
  uncountedBoxData <- reactive({
    uncount(selectedSpecies(), counts)
  })
  
  # make a boxplot
  output$boxPlot <- renderPlot({
    ggplot(uncountedBoxData(), aes(x=species, y=get(input$SpatVars))) +
      geom_boxplot() +
      theme_minimal() +
      ylab(input$SpatVars)
  })
  
  
  
  
  # Data download --------------------------------------------------------------
  
  # Niche Area
  output$downloadNicheAreas <- downloadHandler(
    filename = function() {
      paste("NicheAreas-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NicheAreas, file, row.names = FALSE)
    }
  )
  
  # Niche Range Overlap
  output$downloadNicheOverlaps <- downloadHandler(
    filename = function() {
      paste("NicheRangeOverlap-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(NicheRangeOverlap, file, row.names = FALSE)
    }
  )
  
  # Geospatial Metadata
  output$downloadGeospatial <- downloadHandler(
    filename = function() {
      paste("GeospatialMetadata-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(GeospatialMetadata, file, row.names = FALSE)
    }
  )
   
  
  
}


# Run the app ------------------------------------------------------------------
shinyApp(ui = ui, server = server)

