# @file create_country_data.R
# @brief Convert WHO-supplied CSV files into format required by particle filter

# generated 2 objects: data -- is a dataframe with input values in the columns,
# rows are years sia.info -- is a matrix with rows equal to the age-classes and
# columns for each year each value is the sia impact in that age-class and
# year.

genInput <- function(country) {

    n.ages <- 100 # Total number of age classes for which to generate SIA data

    cases <- read.csv("input-cases.csv", row.names = 1)
    # Columns for input-cases.csv are ISO_code, followed by years, starting
    # with 1980.  We use the number of columns to infer the number of years for
    # which we have data.
    n.yrs <- dim(cases)[2]

    births <- read.csv("input-births.csv", row.names = 1)
    # Columns for input-births.csv are ISO3_code, followed by years, starting
    # with 1980.

    pop <- read.csv("input-populations.csv", row.names = 1)
    # Columns for input-populations.csv begin with an unlabled ISO code column,
    # followed by years, starting with 1980.

    mcv1 <- read.csv("input-mcv1.csv", row.names = 1)
    # Columns for input-mcv1.csv are identical to input-births.csv
    mcv2 <- read.csv("input-mcv2.csv", row.names = 1)
    # Columns for input-mcv2.csv are identical to input-births.csv

    sia.data <- read.csv("input-sia.csv")
    # SIAs are reported alphabetically by country, and in ascending order by
    # year
    high.year <- read.csv("input-outbreak_list.csv")
    # input-outbreak_list.csv reports outbreaks by countries and years

    data <- data.frame(
                       C = unlist(cases[which(rownames(cases) == country), ]),
                       N = unlist(pop[which(rownames(pop) == country), ]), 
                       B = unlist(births[which(rownames(births) == country), ]),
                       mcv1.true = unlist(mcv1[which(rownames(mcv1) == country), ]), 
                       mcv2.true = unlist(mcv2[which(rownames(mcv2) == country), ]),
                       high.year = rep(1, n.yrs)
    )

    # Create SIA information matrix 
    yrs <- dim(data)[1]
    sia.facts = sia.data[which(sia.data$iso == country), ]
    years = c(sia.facts$year) - 1979

    sia.info = matrix(0, nrow = n.ages, ncol = yrs)
    if (length(years) > 0) {

        for (y in 1:length(years)) {
            agegroups = seq(ceiling(sia.facts$agemin[y]), ceiling(sia.facts$agemax[y]))
            sia.info[agegroups, years[y]] = 1

            if(  (ceiling(sia.facts$agemin[y]) == 1 ) || (ceiling(sia.facts$agemin[y]) == 0 )  ){ 
             sia.info[1, years[y]] = (1 - sia.facts$agemin[y]) * 0.85 * sia.facts$coverage[y]/100 
            }


            sia.info[max(agegroups), years[y]] = sia.facts$agemax[y]/ceiling(sia.facts$agemax[y])
            

            sia.info[-1, years[y]] = sia.info[-1, years[y]] * 0.93 * sia.facts$coverage[y]/100 
        }
    }

    data$SIA <- apply(sia.info, 2, max)

    # draw high years from list
    high <- rep(1, dim(data)[1])
    high[high.year$year[which(high.year$iso == country)] - 1979] <- 0
    data$high.years <- high
    data$high.year <- NULL

    write.table(data, file = paste0(country, ".txt"), row.names = FALSE)
    data <- transform(data, SIA = c(0, SIA[-nrow(data)]))
    write.table(data, file = paste0(country, "shift.txt"), row.names = FALSE)
    save(sia.info, file = paste0("sia", country, ".info"))

    write.table(sia.info, file = paste0("sia", country, ".txt"), row.names = TRUE, sep = ",")

    sia.info <- sia.info[, -n.yrs]
    empty <- rep(0, 25)
    new_m <- cbind(c(empty), sia.info)
    write.table(new_m, file = paste0("sia", country, "shift.txt"), row.names = TRUE, sep = ",")
}


# Output files will be generated for the following list of countries:
countries <- c("AFG", "AGO", "ALB", "AND", "ARE", "ARG", "ARM", "ATG", "AUS", "AUT",
    "AZE", "BDI", "BEL", "BEN", "BFA", "BGD", "BGR", "BHR", "BHS", "BIH", "BLR",
    "BLZ", "BOL", "BRA", "BRB", "BRN", "BTN", "BWA", "CAF", "CAN", "CHE", "CHL",
    "CHN", "CIV", "CMR", "COD", "COG", "COK", "COL", "COM", "CPV", "CRI", "CUB",
    "CYP", "CZE", "DEU", "DJI", "DMA", "DNK", "DOM", "DZA", "ECU", "EGY", "ERI",
    "ESP", "EST", "ETH", "FIN", "FJI", "FRA", "FSM", "GAB", "GBR", "GEO", "GHA",
    "GIN", "GMB", "GNB", "GNQ", "GRC", "GRD", "GTM", "GUY", "HND", "HRV", "HTI",
    "HUN", "IDN", "IND", "IRL", "IRN", "IRQ", "ISL", "ISR", "ITA", "JAM", "JOR",
    "JPN", "KAZ", "KEN", "KGZ", "KHM", "KIR", "KNA", "KOR", "KWT", "LAO", "LBN",
    "LBR", "LBY", "LCA", "LKA", "LSO", "LTU", "LUX", "LVA", "MAR", "MCO", "MDA",
    "MDG", "MDV", "MEX", "MHL", "MKD", "MLI", "MLT", "MMR", "MNE", "MNG", "MOZ",
    "MRT", "MUS", "MWI", "MYS", "NAM", "NER", "NGA", "NIC", "NIU", "NLD", "NOR",
    "NPL", "NRU", "NZL", "OMN", "PAK", "PAN", "PER", "PHL", "PLW", "PNG", "POL",
    "PRK", "PRT", "PRY", "QAT", "ROU", "RUS", "RWA", "SAU", "SDN", "SEN", "SGP",
    "SLB", "SLE", "SLV", "SMR", "SOM", "SRB", "STP", "SUR", "SVK", "SVN", "SWE",
    "SWZ", "SYC", "SYR", "TCD", "TGO", "THA", "TJK", "TKM", "TLS", "TON", "TTO",
    "TUN", "TUR", "TUV", "TZA", "UGA", "UKR", "URY", "USA", "UZB", "VCT", "VEN",
    "VNM", "VUT", "WSM", "YEM", "ZAF", "ZMB", "ZWE")


for (i in 1:length(countries)) {
    print(countries[i])
    genInput(countries[i])
}
