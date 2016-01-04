\name{forestTraits}
\docType{data}
\alias{forestTraits}
\title{Trait values for Forest Inventory & Analysis Regions}
\description{Weighted averages by species biomass.  The community weighted mean
values are continuous abundance \code{'CA'} for seed weight, wood density, maximum
height, leaf nitrogen, leaf phosphorus, specific leaf area (SLA), ordinal classes
\code{'OC'} for shade tolerance, drougt tolerance, and flood tolerance, fractional
composition \code{'FC'} for three leaf types (broadleaf deciduous, broadleaf 
evergreen, and needleleaf evergreen).  Dioecy and ringporous are binary traits and,
 thus composition data at the plot scale, but close to zero. They are thus modeled
 as continuous and censored at 0 and 1.}
\usage{forestTraits}
\format{
\tabular{lll}{
  xdata \tab data.frame \tab columns for design matrix\cr
  xdata1 \tab temp \tab winter temperature (standardized)\cr
  xdata2 \tab deficit \tab degree hrs at P < PET\cr
  xdata3 \tab moisture \tab posterior est from FIA index\cr
  xdata4 \tab u1 \tab sin(slope) \cr
  xdata5 \tab u2, u3 \tab u1*sin(aspect), u1*cos(aspect)\cr
  xdata6 \tab stdage \tab FIA stand age\cr
  xdata7 \tab soilFactor \tab soil types are factors\cr
  y \tab traits \tab community-weighted means\cr
  y1 \tab gmPerSeed, gmPerCm, maxHt \tab continuous, standardized\cr
  y2 \tab leafN, leafP, SLA \tab continuous, standardized\cr
  y3 \tab shade, drought, flood \tab ordinal classes 1 - 5\cr
  y4 \tab broaddeciduous, ... \tab three leaf types, composition data\cr
  y5 \tab dioecious, ringPorous \tab binary, two classes each\cr
  typeNames \tab character vector \tab one column per trait\cr
  }
}
\source{Clark, J.S. (in review). Why species tell more about traits than vice versa: 
Predictive analysis to connect species, traits, and environment.}
\keyword{datasets}
