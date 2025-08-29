\name{fingerPro-package}
\alias{fingerPro-package}
\alias{fingerPro}
\docType{package}
\title{
  Sediment Source Fingerprinting
}
\description{
Soil erosion poses significant challenges for food production and leads to reservoir siltation globally. Effective control strategies require information on sediment, nutrient, and pollutant transport. Traditional monitoring techniques often struggle to provide accurate source estimates, but sediment source fingerprinting has proven to be a valuable tool. This package facilitates sediment provenance quantification in catchments, essential for developing erosion prevention and management plans.

The package implements a mixing model algorithm to apportion sediment sources from composite samples. Functions are provided to:
  \itemize{
    \item Characterize and preprocess data, selecting optimal tracers and unraveling multiple solutions.
    \item Unmix sediment samples to estimate source contributions.
    \item Assess source variability effects.
    \item Visualize and export results.
  }
}
\author{
Ivan Lizaga, Borja Latorre, Leticia Gaspar, Ana Navas

Maintainer: Ivan Lizaga <lizaga.ivan10@gmail.com>
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{ Sediment source fingerprinting package }
\seealso{
  GitHub repository:
  \url{https://github.com/eead-csic-eesa/fingerPro}
  
  VIDEO TUTORIALS:
  \url{https://www.youtube.com/watch?v=LcrM_vLOa_I&ab_channel=UnmixingScience}
}
\examples{
# Created by Ivan Lizaga on 26/06/2024

# If you want to use your own data,
# I recommend following the format of the datasets included in the package
# setwd("the directory that contains your dataset")
# data <- read.table('your dataset.csv', header = TRUE, sep = ',')

# Example using the data included in the fingerPro package
# Load the dataset "catchment"
install.packages("fingerPro")
library(fingerPro)

# "catchment": Dataset from a Mediterranean catchment
# for demonstrating package functions. Contains high-quality radionuclides
# and geochemistry data, representing four sources (AG, PI, PI1, SS) and one mixture.

# Load the 'catchment' dataset for testing graphical functions
data <- catchment

# Display boxplots and a correlation graph of the dataset
boxPlot(data, tracers = 1:6, ncol = 3)
correlationPlot(data, columns = 1:7, mixtures = TRUE)

# Check how tracers discriminate between potential sources
 LDAPlot(data[, 1:10], P3D = FALSE, text = FALSE)
 LDAPlot(data[, 1:10], P3D = FALSE, text = TRUE)  # Include point information
 LDAPlot(data[, 1:10], P3D = TRUE, text = FALSE)   # 3D LDA plot

# Interactive 3D LDA plot
# LDAPlot(data[, 1:10], P3D = TRUE, text = TRUE, interactive = TRUE,
#       colors = c("#993300", "#33CC00", "#336600", "#9933CC"))

# Display a PCA graph
PCAPlot(data, colors = c("#993300", "#33CC00", "#336600", "#9933CC", "blue"))

# Combine PI and PI1 sources based on previous analysis
data$Land_Use[data$Land_Use == 'PI1'] <- 'PI'

# Check discriminant capacity of the first 10 tracers again
LDAPlot(data[, 1:10], P3D = FALSE, text = TRUE, colors = c("#993300", "#336600", "#9933CC"))
}
