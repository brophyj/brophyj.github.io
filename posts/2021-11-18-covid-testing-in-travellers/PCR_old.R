library(rvest)
library(janitor)
library(tidyverse)


content <- read_html("https://www.canada.ca/en/public-health/services/diseases/coronavirus-disease-covid-19/testing-screening-contact-tracing/summary-data-travellers.htm")
tables <- content %>% html_table(fill = TRUE)
vacc_august <- tables[[6]][2:9,1:5]  %>% clean_names()
vacc_august$air_2 <- as.numeric(gsub("%","",vacc_august$air_2))
vacc_august$land_2 <- as.numeric(gsub("%","",vacc_august$land_2))
vacc_august$air <- as.numeric(gsub(",","",vacc_august$air))
vacc_august$land <- as.numeric(gsub(",","",vacc_august$land))
names(vacc_august) <- c("Date", "air", "air_p", "land", "land_p")
vacc_august$Date <- gsub(",.*","",vacc_august$Date)

vacc_july <- tables[[5]][2:6,1:5]  %>% clean_names()
vacc_july$air_2 <- as.numeric(gsub("%","",vacc_july$air_2))
vacc_july$land_2 <- as.numeric(gsub("%","",vacc_july$land_2))
vacc_july$air <- as.numeric(gsub(",","",vacc_july$air))
vacc_july$land <- as.numeric(gsub(",","",vacc_july$land))
names(vacc_july) <- c("Date", "air", "air_p", "land", "land_p")
vacc_july$Date <- gsub(",.*","",vacc_july$Date)

vacc <- rbind(vacc_july, vacc_august) %>%
  mutate(air_cases = round(air_p*air/100,0), land_cases = round(land_p*land/100,0), type = "vaccinated")
vair_avg <- sum(vacc$air_p)/sum(vacc$air)


unvacc_august <- tables[[4]][2:9,1:5]  %>% clean_names()
unvacc_august$air_2 <- as.numeric(gsub("%","",unvacc_august$air_2))
unvacc_august$land_2 <- as.numeric(gsub("%","",unvacc_august$land_2))
unvacc_august$air <- as.numeric(gsub(",","",unvacc_august$air))
unvacc_august$land <- as.numeric(gsub(",","",unvacc_august$land))
names(unvacc_august) <- c("Date", "air", "air_p", "land", "land_p")
unvacc_august$Date <- gsub(",.*","",unvacc_august$Date)

unvacc_july <- tables[[3]][2:6,1:5]  %>% clean_names()
unvacc_july$air_2 <- as.numeric(gsub("%","",unvacc_july$air_2))
unvacc_july$land_2 <- as.numeric(gsub("%","",unvacc_july$land_2))
unvacc_july$air <- as.numeric(gsub(",","",unvacc_july$air))
unvacc_july$land <- as.numeric(gsub(",","",unvacc_july$land))
names(unvacc_july) <- c("Date", "air", "air_p", "land", "land_p")
unvacc_july$Date <- gsub(",.*","",unvacc_july$Date)

unvacc <- rbind(unvacc_july, unvacc_august) %>%
  mutate(air_cases = round(air_p*air/100,0), land_cases = round(land_p*land/100,0), type = "unvaccinated")
uvair_avg <- sum(unvacc$air_p)/sum(unvacc$air)

dat <- rbind(vacc, unvacc)
dat$Date<- gsub("-.*","",dat$Date)
dat$Date<- gsub(" ","/",dat$Date)
dat$Date <- as.Date(dat$Date, "%B/%d")
dat$comb <- dat$air + dat$land
dat <- dat %>%
  mutate(comb_cases = land_cases + air_cases, comb_p = comb_cases / comb)

early <- tables[[2]][3:6,1:5]  %>% clean_names()
early$air_2 <- as.numeric(gsub("%","",early$air_2))
early$land_2 <- as.numeric(gsub("%","",early$land_2))
early$air <- as.numeric(gsub(",","",early$air))
early$land <- as.numeric(gsub(",","",early$land))
names(early) <- c("Date", "air", "air_p", "land", "land_p")
early$Date <- gsub(",.*","",early$Date)

g <- ggplot(dat, aes(x=Date, y=comb_p, group=type)) +
  geom_line(aes(linetype=type, color=type)) +
  geom_point()  + 
  labs(title="Percentage of + PCR tests among Canadian travellers",
       x ="Week", y = "Percentage  postive tests", 
       caption = "https://www.canada.ca/en/public-health/services\n/diseases/coronavirus-disease-covid-19\n/testing-screening-contact-tracing/\nsummary-data-travellers.htm") +
  geom_smooth() +
  scale_y_continuous(name="Percentage  postive tests", labels = function(x) paste0(x, "%")) +
  theme_bw()

g
pdf("PCR.pdf")
g
dev.off()

png("PCR.png")
g
dev.off()


(sum(dat[1:13,5]) + sum(dat[1:13,6])) / (sum(dat[1:13,2]) + sum(dat[1:13,4]))
(sum(dat[1:13,2]) + sum(dat[1:13,4])) * 200
sum(dat[1:13,6]) + sum(dat[1:13,7])
sum(dat[1:13,10]) / sum(dat[1:13,9])


