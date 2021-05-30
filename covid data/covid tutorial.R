# COVID-19 Tutorial
#------------------------
library(tmap)
# IMPORT RAW DATA: Johns Hopkins Github data
confirmedraw <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
str(confirmedraw) # Check latest date at the end of data
deathsraw <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
recoveredraw <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")
# Note differences in the number of rows/columns


# DATA CLEANING: To create country level and global combined data
# Convert each data set from wide to long AND aggregate at country level
library(tidyr)
library(dplyr)
confirmed <- confirmedraw %>% gather(key="date", value="confirmed", -c(Country.Region, Province.State, Lat, Long)) %>% group_by(Country.Region, date) %>% summarize(confirmed=sum(confirmed))
deaths <- deathsraw %>% gather(key="date", value="deaths", -c(Country.Region, Province.State, Lat, Long)) %>% group_by(Country.Region, date) %>% summarize(deaths=sum(deaths))
recovered <- recoveredraw %>% gather(key="date", value="recovered", -c(Country.Region, Province.State, Lat, Long)) %>% group_by(Country.Region, date) %>% summarize(recovered=sum(recovered))
summary(confirmed)

# Final data: combine all three
country <- full_join(confirmed, deaths) %>% full_join(recovered)

# Date variable
# Fix date variable and convert from character to date
str(country) # check date character
country$date <- country$date %>% sub("X", "", .) %>% as.Date("%m.%d.%y")
str(country) # check date Date
# Create new variable: number of days
country <- country %>% group_by(Country.Region) %>% mutate(cumconfirmed=cumsum(confirmed), days = date - first(date) + 1)

# Aggregate at world level
world <- country %>% group_by(date) %>% summarize(confirmed=sum(confirmed), cumconfirmed=sum(cumconfirmed), deaths=sum(deaths), recovered=sum(recovered)) %>% mutate(days = date - first(date) + 1)
# Extract specific country: Italy
italy <- country %>% filter(Country.Region=="Italy")



# SUMMARY STATISTICS
summary(country)
by(country$confirmed, country$Country.Region, summary)
by(country$cumconfirmed, country$Country.Region, summary)
by(country$deaths, country$Country.Region, summary)
by(country$recovered, country$Country.Region, summary)
summary(world)
summary(italy)



# GRAPHS
# Barchart of cases over time
library(ggplot2)
# World confirmed
ggplot(world, aes(x=date, y=confirmed)) + geom_bar(stat="identity", width=0.1) +
  theme_classic() +
  labs(title = "Covid-19 Global Confirmed Cases", x= "Date", y= "Daily confirmed cases") +
  theme(plot.title = element_text(hjust = 0.5))

# Italy confirmed
ggplot(italy, aes(x=date, y=confirmed)) + geom_bar(stat="identity", width=0.1) +
  theme_classic() +
  labs(title = "Covid-19 Confirmed Cases in Italy", x= "Date", y= "Daily confirmed cases") +
  theme(plot.title = element_text(hjust = 0.5))

# World confirmed, deaths and recovered
str(world)
world %>% gather("Type", "Cases", -c(date, days)) %>%
  ggplot(aes(x=date, y=Cases, colour=Type)) + geom_bar(stat="identity", width=0.2, fill="white") +
  theme_classic() +
  labs(title = "Covid-19 Global Cases", x= "Date", y= "Daily cases") +
  theme(plot.title = element_text(hjust = 0.5))


# Line graph of cases over time
# World confirmed
ggplot(world, aes(x=days, y=confirmed)) + geom_line() +
  theme_classic() +
  labs(title = "Covid-19 Global Confirmed Cases", x= "Days", y= "Daily confirmed cases") +
  theme(plot.title = element_text(hjust = 0.5))
# Ignore warning

# World confirmed with counts in log10 scale
ggplot(world, aes(x=days, y=confirmed)) + geom_line() +
  theme_classic() +
  labs(title = "Covid-19 Global Confirmed Cases", x= "Days", y= "Daily confirmed cases  (log scale)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(trans="log10")

# World confirmed, deaths and recovered
str(world)
world %>% gather("Type", "Cases", -c(date, days)) %>%
  ggplot(aes(x=days, y=Cases, colour=Type)) + geom_line() +
  theme_classic() +
  labs(title = "Covid-19 Global Cases", x= "Days", y= "Daily cases") +
  theme(plot.title = element_text(hjust = 0.5))

# Confirmed by country for select countries with counts in log10 scale
countryselection <- country %>% filter(Country.Region==c("US", "Italy", "China", "France", "United Kingdom", "Germany"))
ggplot(countryselection, aes(x=days, y=confirmed, colour=Country.Region)) + geom_line(size=1) +
  theme_classic() +
  labs(title = "Covid-19 Confirmed Cases by Country", x= "Days", y= "Daily confirmed cases (log scale)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(trans="log10")

# Matrix of line graphs of confirmed, deaths and recovered for select countries in log10 scale
str(countryselection)
countryselection %>% gather("Type", "Cases", -c(date, days, Country.Region)) %>%
  ggplot(aes(x=days, y=Cases, colour=Country.Region)) + geom_line(size=1) +
  theme_classic() +
  labs(title = "Covid-19 Cases by Country", x= "Days", y= "Daily cases (log scale)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(trans="log10") +
  facet_grid(rows=vars(Type))

#---------------------------------------

## Map
countrytotal <- country %>% group_by(Country.Region) %>% summarize(cumconfirmed=sum(confirmed), cumdeaths=sum(deaths), cumrecovered=sum(recovered))

# Basemap from package tmap
library(tmap)
data(World)
class(World)

# Combine basemap data to covid data
countrytotal$Country.Region[!countrytotal$Country.Region %in% World$name]
list <- which(!countrytotal$Country.Region %in% World$name)
countrytotal$country <- as.character(countrytotal$Country.Region)
countrytotal$country[list] <-
  c("Andorra", "Antigua and Barbuda", "Bahrain",
    "Barbados", "Bosnia and Herz.", "Myanmar",
    "Cape Verde", "Central African Rep.", "Congo",
    "Dem. Rep. Congo", "Czech Rep.", "Diamond Princess",
    "Dominica", "Dominican Rep.", "Eq. Guinea",
    "Swaziland", "Grenada", "Holy See",
    "Korea", "Lao PDR", "Liechtenstein",
    "Maldives", "Malta", "Mauritius",
    "Monaco", "MS Zaandam", "Macedonia",
    "Saint Kitts and Nevis", "Saint Lucia", "Saint Vincent and the Grenadines",
    "San Marino", "Sao Tome and Principe", "Seychelles",
    "Singapore", "S. Sudan", "Taiwan",
    "United States", "Palestine", "W. Sahara")
countrytotal$Country.Region[!countrytotal$country %in% World$name]
World$country <- World$name
worldmap <- left_join(World, countrytotal, by="country")
worldmap$cumconfirmed[is.na(worldmap$cumconfirmed)] <- 0

# Map
ggplot(data = worldmap) + geom_sf(aes(fill=cumconfirmed), color="black") +
  ggtitle("World Map of Confirmed Covid Cases",
          subtitle="Total Cases on April 20, 2020") +
  theme_bw()

#--------------------------------------------

## MODEL

# DATA FOR REGION: Ontario
ontario <- confirmedraw %>% filter(Province.State=="Ontario") %>% gather(key="date", value="confirmed", -c(Country.Region, Province.State, Lat, Long)) %>%  mutate(cumconfirmed=cumsum(confirmed))
ontario$date <- ontario$date %>% sub("X", "", .) %>% as.Date("%m.%d.%y")

# SIR FUNCTION
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S/N
    dI <- beta * I * S/N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}

# CREATE A VECTOR OF DAILY CUMULATIVE INCIDENCE NUMBERS OF CANADA FROM START DATE
infected <- ontario %>% filter(confirmed>0) %>% pull(cumconfirmed)

# Create an incrementing Day vector the same length as our cases vector
day <- 1:(length(infected))
N <- 14446515
# now specify initial values for S, I and R
init <- c(S = N - infected[1], I = infected[1], R = 0)

RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = day, func = SIR, parms = parameters)
  fit <- out[, 3]
  sum((infected - fit)^2)
}

# now find the values of beta and gamma that give the
# smallest RSS, which represents the best fit to the data.
# Start with values of 0.5 for each, and constrain them to
# the interval 0 to 1.0
optimization <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", lower = c(0,0), upper = c(1, 1))

# check for convergence
optimization$message

# Optimization Parameters
opt_par <- setNames(optimization$par, c("beta", "gamma"))
opt_par

# Reproduction Number
R0 <- opt_par[1]/opt_par[2]
R0

# PREDICTION
# time in days for predictions
t <- 1:150
# get the fitted values from our SIR model
fittedcum <- data.frame(ode(y = init, times = t, func = SIR, parms = opt_par))
# add a Date column and join the observed incidence data
fittedcum <- fittedcum %>%
  mutate(date = as.Date(startdate) + t - 1, province = "Ontario") %>%
  left_join(ontario %>% select(date, cumconfirmed))

# plot the data
ggplot(fittedcum, aes(x = date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_point(aes(y = cumconfirmed), colour = "orange")
labs(y = "Cumulative incidence", x="Date",
     title = "COVID-19 fitted vs observed cumulative incidence, Ontario",
     subtitle = "(red=fitted incidence from SIR model, orange=observed incidence)")

# plot the data
ggplot(fittedcum, aes(x = date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_line(aes(y = S), colour = "black") +
  geom_line(aes(y = R), colour = "green") +
  geom_point(aes(y = cumconfirmed), colour = "orange") +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Persons", title = "COVID-19 fitted vs observed cumulative incidence, Ontario province") +
  scale_colour_manual(name = "",
                      values = c(red = "red", black = "black", green = "green", orange = "orange"),
                      labels = c("Susceptible", "Recovered", "Observed incidence", "Infectious")) +
  scale_y_continuous(trans="log10")
