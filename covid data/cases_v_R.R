cases_infections <- read.csv("/Users/brycemorsky/Desktop/New projects/Norms and pandemics/code/norms and pandemics/covid data/us_daily_cases_and_r.csv")

infections = cases_infections[1:400,3];
cases = cases_infections[1:400,4];

alabama = cases_infections[1:400,3:4];

library(ggplot2)
library(dplyr)

filter(!is.na(cases_infections[420:500,:]))

ggplot(cases_infections[1:400,3:4], aes(x=New.Cases, y=Estimated.Effective.R)) + geom_point()

cor(infections,cases,method = c("pearson", "kendall", "spearman"))