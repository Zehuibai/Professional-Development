# Install and import required libraries
require(shiny)
require(ggplot2)
require(leaflet)
require(tidyverse)
require(httr)
require(scales)
# Import model_prediction R which contains methods to call OpenWeather API
# and make predictions
source("model_prediction.R")


test_weather_data_generation<-function(){
  #Test generate_city_weather_bike_data() function
  city_weather_bike_df<-generate_city_weather_bike_data()
  stopifnot(length(city_weather_bike_df)>0)
  print(head(city_weather_bike_df))
  return(city_weather_bike_df)
}

# Create a RShiny server
shinyServer(function(input, output){
  # Define a city list
  city <- c("Seoul", "Suzhou", "New York", "Paris", "London")
  # Define color factor
  color_levels <- colorFactor(c("green", "yellow", "red"), 
                              levels = c("small", "medium", "large"))
  city_weather_bike_df <- test_weather_data_generation()
  
  # Create another data frame called `cities_max_bike` with each row 
  # contains city location info and max bike prediction for the city
  cities <- city_weather_bike_df %>% group_by(CITY_ASCII)
  cities_max_bike <- cities %>% filter(BIKE_PREDICTION == max(BIKE_PREDICTION))
  
  output$city_bike_map <- 
    renderLeaflet({    
      leaflet() %>%
      addTiles() %>%
      addCircleMarkers(data = cities_max_bike, lng = cities_max_bike$LNG, lat = cities_max_bike$LAT, 
                      popup = cities_max_bike$LABEL,
                      color = ~color_levels(cities_max_bike$BIKE_PREDICTION_LEVEL), 
                      radius = ~ifelse(cities_max_bike$BIKE_PREDICTION_LEVEL=='small', 6, 12))
      # Complete this function to render a leaflet map
    })
  
  
    observeEvent(input$city_dropdown, {
    if(input$city_dropdown != 'All') {
      
      leafletProxy("city_bike_map") %>% clearShapes()
      index = which(cities_max_bike$CITY_ASCII == input$city_dropdown)
      leafletProxy("city_bike_map") %>%
      addCircles(lng = cities_max_bike$LNG[index], 
                lat = cities_max_bike$LAT[index],
                popup=cities_max_bike$DETAILED_LABEL[index])
                
      new_cities <- city_weather_bike_df %>% filter(city_weather_bike_df$CITY_ASCII == input$city_dropdown)
      
      output$temp_line <- renderPlot({
        ggplot(new_cities, aes(x = FORECASTDATETIME, y=TEMPERATURE)) +
        geom_text(label=new_cities$TEMPERATURE) + 
        geom_point()+
        geom_line(color='yellow', group = 1)
      })
      
      output$bike_line <- renderPlot({
        ggplot(new_cities, aes(x = new_cities$FORECASTDATETIME, y=new_cities$BIKE_PREDICTION)) +
        geom_text(label=new_cities$BIKE_PREDICTION) +
        geom_point()+
        geom_line(color='red',group = 1)
      })
      
      output$bike_data_output <- renderText({
        paste("Time=", input$plot_clicks$x, "\nBikeCountPred=", input$plot_click$y)
      })
      
      output$humidity_pred_chart <- renderPlot({
        ggplot(new_cities, aes(x = new_cities$HUMIDITY, y = new_cities$BIKE_PREDICTION)) +
        geom_point()+
        geom_smooth(formula = y ~ poly(x,4), method = "lm")
      })
                      
    }
      else {
        leafletProxy("city_bike_map")}
      })
  # Observe drop-down event
  
  # Then render output plots with an id defined in ui.R
  
  # If All was selected from dropdown, then render a leaflet map with circle markers
  # and popup weather LABEL for all five cities
  
  # If just one specific city was selected, then render a leaflet map with one marker
  # on the map and a popup with DETAILED_LABEL displayed
  
})
