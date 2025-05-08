#### Function for reading in and formatting NLDAS data
### Author: Jacob Schaperow 
## Data accessed from "https://disc.gsfc.nasa.gov/information/tools?title=Hydrology%20Data%20Rods"

get_NLDASv20_datarod <- function(start_date, end_date, lat, lon, var)
{
  base_url <- "https://hydro1.gesdisc.eosdis.nasa.gov/daac-bin/access/timeseries.cgi"
  
  full_url <- paste0(base_url,
                     "?variable=NLDAS2:NLDAS_FORA0125_H_v2.0:", var,
                     "&startDate=", start_date, "T00&",
                     "endDate=", end_date, "T23&",
                     "location=GEOM:POINT(", lon, ",%20", lat, ")",
                     "&type=asc2"
  )
  print(full_url)
  
  x <- GET(full_url)
  
  # Parse response
  z <- content(x, as = "text")
  z1 <- strsplit(z, "\n")[[1]] # separate data by line
  
  # Get dates and values. First 13 lines are metadata.
  z1 <- z1[14:length(z1)]
  
  # Split each entry by the tab separator into datetime and value
  split_entries <- strsplit(z1, "\t")
  
  # Extract datetime and value columns
  thedatetime <- sapply(split_entries, function(x) x[1]) # First part is datetime
  thevalue <- sapply(split_entries, function(x) as.numeric(x[2])) # Second part is value (convert to numeric)
  
  # Create the data frame
  df <- data.frame(datetime = as.POSIXct(thedatetime, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"), value = thevalue)
  
  return(df)
}