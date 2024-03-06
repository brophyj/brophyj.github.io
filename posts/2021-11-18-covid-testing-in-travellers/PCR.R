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

#new vaccinated
vac <- tables[[4]][2:10,1:5] %>% clean_names()
vac$air <- as.numeric(gsub(",","",vac$air))
vac$land <- as.numeric(gsub(",","",vac$land ))
vac$air_2 <- as.numeric(gsub("%","",vac$air_2))
vac$land_2 <- as.numeric(gsub("%","",vac$land_2))
names(vac) <- c("Date", "air", "air_p", "land", "land_p")
vac$Date <- gsub(",.*","",vac$Date)
vac$Date <- gsub(".*\\(","",vac$Date)
vac <- vac %>%
  mutate(air_n=round(air*air_p/100,0), land_n=round(land*land_p/100,0), comb=air_n+land_n, comb_p= 100*comb/(air+land), type="vac")
vac

#new unvaccinated

unvac <- tables[[3]][2:10,1:5] %>% clean_names()
unvac$air <- as.numeric(gsub(",","",unvac$air))
unvac$land <- as.numeric(gsub(",","",unvac$land ))
unvac$air_2 <- as.numeric(gsub("%","",unvac$air_2))
unvac$land_2 <- as.numeric(gsub("%","",unvac$land_2))
names(unvac) <- c("Date", "air", "air_p", "land", "land_p")
unvac$Date <- gsub(",.*","",unvac$Date)
unvac$Date <- gsub(".*\\(","",unvac$Date)
unvac <- unvac %>%
  mutate(air_n=round(air*air_p/100,0), land_n=round(land*land_p/100,0), comb=air_n+land_n, comb_p= 100*comb/(air+land), type="unvac")
unvac

dat <- rbind(vac, unvac)

dat$Date <- factor(dat$Date, levels=unique(dat$Date)) # needed to stop ggplot sorting alphabetically

g <- ggplot(dat, aes(x=Date, y=comb_p, group=type)) +
  geom_line(aes(linetype=type, color=type)) +
  geom_point()  + 
  labs(title="Percentage of + PCR tests among Canadian travellers",
       x ="Time period", y = "Percentage  postive tests", 
       caption = "https://www.canada.ca/en/public-health/services\n/diseases/coronavirus-disease-covid-19\n/testing-screening-contact-tracing/\nsummary-data-travellers.htm") +
  geom_smooth() +
  scale_y_continuous(name="Percentage  postive tests", labels = function(x) paste0(x, "%")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1))

g
pdf("PCR.pdf")
g
dev.off()

png("PCR.png")
g
dev.off()


paste("Total number vaccinated air travellers tested is ", sum(dat[1:9,2]), " and of these ", sum(dat[1:9,6]), " or ", round(mean(dat$air_p[1:9]),2), "% tested positive")
paste("Total number vaccinated land travellers tested is ", sum(dat[1:9,4]), " and of these ", sum(dat[1:9,7]), " or ", round(mean(dat$land_p[1:9]),2), "% tested positive")
paste("Total number vaccinated travellers tested is ", sum(dat[1:9,2]) + sum(dat[1:9,4]), " and of these ", sum(dat[1:9,8]), " or ", round(mean(dat$comb_p[1:9]),2), "% tested positive")

