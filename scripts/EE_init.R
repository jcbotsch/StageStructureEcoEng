#====load packages=====
library(tidyverse)
library(lubridate)
library(lme4)
library(car)

theme_set(theme_classic())
#Define function to give NA if all NAs else give mean
meanna <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    mean(x, na.rm = TRUE)
  }
}

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


resp_avg <- inc_full %>% 
  filter(sampledate>as.Date("2017-05-24"),
         dark_light == 0) %>% 
  group_by(sampledate) %>% 
  summarise(do_gm2hr = mean(do_gm2hr),
            wtemp = mean(wtemp)) %>% 
  mutate(light50 = 0,
         light_logger = 0) 

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


resp.data.avg <- resp_avg %>% 
  left_join(midges %>% 
              group_by(sampledate, instar) %>% 
              summarise(larvae = meanna(larvae),
                        ninstar = meanna(ninstar))) %>% 
  filter(!is.na(larvae))

data <- nep %>% 
  left_join(midges) %>% 
  filter(!is.na(larvae)) %>% 
  full_join(resp.data.avg)




data %>% 
  ggplot(aes(x = light50, y = do_gm2hr, size = ninstar, col = year(sampledate) + yday(sampledate)/365))+
  facet_wrap(~instar)+
  geom_point(alpha = 0.5)+
  labs(color = "Date",
       size = "Larvae")+
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

mod2.0 <- lmer(do_gm2hr~log1p(larvae)+(1|sampledate), data = data %>% mutate(sampledate = as.character(sampledate), year= as.character(year(sampledate)))) 

summary(mod2.0)
Anova(mod2.0, test.statistic = "F")

data %>% 
  group_by(sampledate) %>%
  mutate(do_var = do_gm2hr-mean(do_gm2hr)) %>% 
  ggplot(aes(x = log1p(larvae), y = do_gm2hr))+
  geom_point()+
  geom_smooth(method = "lm")

data %>% 
  group_by(sampledate) %>%
  mutate(do_var = do_gm2hr-mean(do_gm2hr),
         larvvar = (larvae-mean(larvae))/mean(larvae)) %>% 
  ggplot(aes(x = larvvar, y = do_var))+
  geom_point()+
  geom_smooth(method = "lm")


data %>% 
  ggplot(aes(x = log1p(ninstar), y = do_gm2hr))+
  facet_wrap(~instar)+
  geom_point()+
  geom_smooth(method = "lm")

datavg <- data %>% 
  mutate(instar = as.numeric(str_remove(instar,"s")),
         star = instar*ninstar) %>% 
  group_by(sampledate, coreid, do_gm2hr, larvae) %>% 
  summarise(avg_instar = sum((ninstar*instar)/(larvae), na.rm = TRUE))

mod2.3 <- lmer(do_gm2hr~avg_instar*larvae+ (1|sampledate), data = datavg) 

summary(mod2.3)

Anova(mod2.3, test.statistic = "F", type = 3)
Anova(mod2.3, test.statistic = "F", type = 2)


mod2 <- lmer(do_gm2hr~s1+s2+s3+s4 + (1|sampledate), data = datawide) 

summary(mod2)

Anova(mod2, test.statistic = "F")


mod2.1 <- lmer(do_gm2hr~log1p(s1)+log1p(s2)+log1p(s3)+log1p(s4) + (1|sampledate), data = datawide) 

summary(mod2.1)

Anova(mod2.1)


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


#====Fit midge P-I model using NLS=====
dz <- data %>% 
  # filter(!is.na(coreid)) %>% 
  select(sampledate, light50, larvae, do_gm2hr) %>% 
  unique() %>% 
  mutate(light.z = (light50+1)/ max(light50+1),
         midge.z = log(larvae+1)/max(log(larvae+1)),
         nep.z = (do_gm2hr)/ max(do_gm2hr))


skimr::skim(dz)


mod <- nls(nep.z~( (me*midge.z+ beta)*light.z)/(k + light.z) + (d + med*midge.z),
    data = dz,
    start = c( beta = 0.01, k = 0.5, me = 0.1, d = -0.1, med = -0.1),
    lower = c(0, 0, 0, -50, -50),
    upper = c(50, 50,50, 0, 50),
    algorithm = "port")

summary(mod)
plot(mod)

(coef(mod)["beta"] + coef(mod)["me"])/coef(mod)["k"]


(coef(mod)["beta"])/coef(mod)["k"]


nd = data.frame(light.z = seq(min(dz$light.z), max(dz$light.z), by = 0.01)) %>% 
  crossing(midge.z = c(min(dz$midge.z), max(dz$midge.z)))
nd$nep <- predict(mod, newdata = nd)

nd %>% 
  ggplot(aes(x = light.z, y = nep,))+
  geom_line(aes( color = factor(midge.z)))+
  geom_point(aes(y = nep.z, fill = midge.z), shape= 21, data = dz)+
  scale_color_viridis_d()+
  scale_fill_viridis_c()



mod.nome <- nls(nep.z~( (beta)*light.z)/(k + light.z) + d,
                data = dz,
                start = c( beta = 0.01,k = 0.5, d = -0.1),
                lower = c(0, 0, -50),
                upper = c(50, 50, 0),
                algorithm = "port")

summary(mod.nome)


nd2 <- data.frame(light.z = seq(min(dz$light.z), max(dz$light.z), by = 0.01), midge.z = "no midge efffect")
nd2$nep <- predict(mod.nome, newdata = nd2)

ndfull <- full_join(nd %>% mutate(midge.z = ifelse(midge.z == 1, "peak midge", "no midge")), nd2)

ndfull %>% 
  ggplot(aes(x = light.z, y = nep, ))+
  geom_point(aes(x = light.z, y = nep.z, fill = midge.z), shape = 21, size = 2, alpha = 0.7, data = dz)+
  geom_line(aes(color = as.character(midge.z)))+
  labs(x = "Scaled Light",
       y = "Scaled NEP",
       color = element_blank(),
       fill = "Scaled Midges")+
  scale_color_viridis_d()+
  scale_fill_viridis_c()

# ggpreview(plot = last_plot(), width = 5, height = 3, units = "in", dpi = 650)


logLik(mod)
logLik(mod.nome)
anova(mod, mod.nome)
AIC(mod, mod.nome)

#====Midge effects on alpha and beta====


mod2 <- nls(nep.z~( (me*midge.z+ beta)*light.z)/(((me*midge.z+ beta)/(alpha+mea*midge.z)) + light.z) + (d + med*midge.z),
           data = dz,
           start = c( beta = 0.1, alpha = 0.1, mea = 0.01, me = 0.1, d = -0.1, med = -0.1),
           lower = c(0, 0, 0, 0, -50, -50),
           upper = c(50, 50,50, 50, 0, 50),
           algorithm = "port")

summary(mod2)
plot(mod2)