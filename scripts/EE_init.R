#====load packages=====
library(tidyverse)
library(lubridate)
library(lme4)
library(car)

theme_set(theme_classic())

#====read data====
inc <- read_csv("data/incubation_core.csv")
inc_raw <- read_csv("data/NEP_21.csv")
inc21 <- inc_raw %>% 
  filter(year(sampledate)==2021)


inc_full <- inc %>% 
  full_join(inc21 %>% 
  filter(dark_light!=-1) %>% 
    mutate(delta_do_rate = (do_final - do_init)/as.numeric(inc_time_end-inc_time_start), #mg O2 L-1 h-1
           do_gm2hr = delta_do_rate * column_depth/100, #converting units to g O2 m-2 hr-1 (1 mg/L = 1g/m3)
           light50 = (light50_init+light50_final)/2,
           wtemp = (wtemp_init+wtemp_final)/2)) %>% 
  select(sampledate, coreid, light_logger, light50, dark_light, wtemp, do_gm2hr)

nep <- inc_full %>% 
  filter(sampledate>as.Date("2017-05-24")) %>% 
  mutate(light50 = ifelse(dark_light == 0, 0, light50),
         light_logger = ifelse(dark_light == 0, 0, light_logger)) %>% 
  select(-dark_light)

unique(nep$sampledate)

rc <- read_csv("data/rc_data.csv")


midges <- rc %>% 
  filter(species_name == "tt") %>% 
  select(sta, sampledate, coreid, larvae, instar, ninstar) %>% 
  mutate(instar = paste0("s", instar)) %>%
  spread(instar, ninstar, fill = 0) %>% 
  select(-sNA) %>% 
  gather(instar, ninstar, s1:s4) %>% 
  filter(sta == 33,
         year(sampledate)>2016)

data <- nep %>% 
  left_join(midges) %>% 
  filter(!is.na(larvae))

data %>% 
  ggplot(aes(x = light50, y = do_gm2hr, color = log1p(ninstar)))+
  facet_wrap(~instar)+
  geom_point()+
  scale_color_viridis_c()


data %>% 
  spread(instar, ninstar, fill = 0) %>% 
  gather(instar, ninstar, s1:s4) %>% 
  ggplot(aes(x = light_logger, y = do_gm2hr, color = log1p(ninstar)))+
  facet_wrap(~instar)+
  geom_point()+
  scale_color_viridis_c()


skimr::skim(data)


datawide <- data %>% 
  spread(instar, ninstar, fill = 0)

skimr::skim(datawide)


mod1 <- lmer(log1p(larvae)~do_gm2hr+(1|sampledate) + (1|year), data = data %>% mutate(sampledate = as.character(sampledate), year= as.character(year(sampledate)))) 

summary(mod1)
Anova(mod1, test.statistic = "F")


mod1.5 <- data %>%
  mutate(sampledate = as.character(sampledate), 
         year= as.character(year(sampledate))) %>% 
  group_by(year, sampledate) %>% 
  summarise(meanlarvae = mean(larvae)) %>%
  group_by(year) %>% 
  mutate(laglarvae = lag(meanlarvae)) %>% 
  full_join(data %>% mutate(sampledate = as.character(sampledate))) %>% 
  lmer(log1p(larvae)~do_gm2hr + log1p(laglarvae) + (1|sampledate), data = .) 

summary(mod1.5)

Anova(mod1.5, test.statistic = "F")

data %>% 
  ggplot(aes(x = do_gm2hr, y = log1p(larvae), color = factor(sampledate)))+
  geom_line(aes(group = sampledate), alpha = 0.2)+
  geom_point()+
  theme(legend.position = "none")


data %>% 
  ggplot(aes(x = do_gm2hr, y = log1p(larvae), color = factor(year(sampledate))))+
  geom_line(aes(group = sampledate), alpha = 0.2)+
  geom_point()+
  labs(color = element_blank())+
  theme(legend.position = "bottom")

mod2 <- lmer(do_gm2hr~log1p(s1)+log1p(s2)+log1p(s3)+log1p(s4) + (1|sampledate), data = datawide) 

summary(mod2)

Anova(mod2)

data %>% 
  ggplot(aes(y = do_gm2hr, x = log1p(ninstar), color = factor(sampledate)))+
  facet_wrap(~instar)+
  geom_point()+
  geom_line(aes(group = sampledate), alpha = 0.2)+
  theme(legend.position = "none")


data %>% 
  ggplot(aes(y = do_gm2hr, x = log1p(ninstar), color = light50))+
  # facet_grid(year(sampledate)~instar)+
  facet_wrap(~instar)+
  geom_point()+
  geom_line(aes(group = sampledate), alpha = 0.2)+
  scale_color_viridis_c()+
  theme(legend.position = "bottom")
