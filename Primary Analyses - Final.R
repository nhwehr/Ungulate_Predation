###### Full analysis for "Interspecific carnivore competition and ungulate predation correlate to predator species richness"

##### Code written by Nathaniel H. Wehr, Hailey M. Boone, Merijn van den Bosch, and Alejandra Zubiria Perez
### Last edited on February 23, 2024 ###

# Data Prep ----------
# Load libraries
library(dplyr)
library(corrplot)
library(grid)
library(glmmTMB)
library(ggplot2)
library(gridExtra)
library(ggeffects)
library(ggpubr)
 
# Load data
data <- read.csv("Raw data/Analyzed Dataset - Ecology - V7.csv", header = T)

## Data manipulation and cleanup
columns <- c("Grizzly_Bear_Present", "Black_Bear_Present", "Mountain_Lion_Present", 
             "Bobcat_Present", "Wolf_Present", "Coyote_Present", "Wolverine_Present")

# Put NA for depredation areas where the predator is not present
for (col in columns) {
  data[[col]][data[[col]] == 0] <- NA 
}

### Invert Dataset for Predator Specific Analyses
## Figure out how many ungulates each predator species has killed and reorganize data
process_predator <- function(data, predator_name, present_column, predation_column) {
  predator <- data %>%
    select(present_column, predation_column, Common_Name, Age_Class, Latitude, Longitude, Monitored, Known_Mortalities, Legal_Harvest, Predation, Number_of_Relevant_Predators, Number_of_Ungulates, Human_Footprint_Index, Percent_Forest_Cover, Average_Terrain_Ruggedness) %>%
    mutate(Predator_Common = predator_name, 
           Predation_Count = !!sym(predation_column)) %>% 
    na.omit() %>% 
    select(-predation_column, -present_column) 
  return(predator)
  }

Grizzly <- process_predator(data, "Grizzly bear", "Grizzly_Bear_Present", "Grizzly_Bear_Predation")
BlackBear <- process_predator(data, "Black bear", "Black_Bear_Present", "Black_Bear_Predation")
MountainLion <- process_predator(data, "Cougar", "Mountain_Lion_Present", "Mountain_Lion_Predation")
Bobcat <- process_predator(data, "Bobcat", "Bobcat_Present", "Bobcat_Predation")
Wolf <- process_predator(data, "Gray wolf", "Wolf_Present", "Wolf_Predation")
Coyote <- process_predator(data, "Coyote", "Coyote_Present", "Coyote_Predation")
Wolverine <- process_predator(data, "Wolverine", "Wolverine_Present", "Wolverine_Predation")

## Combine and Clean
# Calculate proportion of mortality each predator made
preddata <- rbind(Grizzly, BlackBear, MountainLion, Bobcat, Wolf, Coyote, Wolverine) %>% 
  mutate(PredPropMort = Predation_Count / Predation)

# Check data distribution
hist(preddata$Latitude)
hist(preddata$Longitude)
hist(preddata$Monitored)
hist(preddata$Known_Mortalities)
hist(preddata$Legal_Harvest)
hist(preddata$Predation)
hist(preddata$Number_of_Relevant_Predators)
hist(preddata$Number_of_Ungulates)
hist(preddata$Human_Footprint_Index)
hist(preddata$Percent_Forest_Cover)
hist(preddata$Average_Terrain_Ruggedness)
hist(preddata$Predation_Count)
hist(preddata$PredPropMort)

# Correlation check
covcortest <- cor(preddata[,c("Number_of_Relevant_Predators","Number_of_Ungulates", "Predation_Count", "PredPropMort", "Legal_Harvest", "Predation", "Human_Footprint_Index", "Percent_Forest_Cover", "Average_Terrain_Ruggedness")], use = 'complete.obs')
corrplot(covcortest)
covcortest

# Model 1 ---------------
# Import data
alldata <- data

# GGplot theme
mytheme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        strip.background = element_rect("white"),
        strip.text = element_text(size = 12))

alldata$random <- paste0("A", as.numeric(factor(alldata$Longitude))) #gives each study a random ID

preddata <- preddata %>% 
  filter(Age_Class != "Mixed") %>% #remove any sites with mixed age class so only adults and neonates remain 
  mutate(Predator_Common = as.character(Predator_Common),
         random = paste0("A", as.numeric(factor(Longitude))))

# Isolate Ungulates
common_names <- unique(preddata$Common_Name)
df_names <- c("Caribou", "Elk", "Moose", "Muledeer", "Whitetail", "Bighorn", "Pronghorn")

split_data <- lapply(common_names, function(name) {
  preddata %>%
    filter(Common_Name == name) %>%
    na.omit() %>%
    mutate_at(vars(Age_Class, Predator_Common, random), as.factor)
})

list2env(setNames(split_data, paste0(df_names, "1")), envir = .GlobalEnv)

## Bighorn ----------
mod1 <- glmmTMB(PredPropMort ~ Age_Class*Predator_Common + Number_of_Relevant_Predators + (1|random),
                betabinomial(link = "logit"),
                weights = Predation,
                data = Bighorn1)
summary(mod1) 
# Analysis was concluded:
#all models have the same number of relevant predators

## Caribou ------------------
mod1 <- glmmTMB(PredPropMort ~ Age_Class*Predator_Common + Number_of_Relevant_Predators + (1|random),
                betabinomial(link = "logit"),
                weights=Predation,
                data=Caribou1)
summary(mod1)
# analysis was concluded:
#model 3 was better fitting for Caribou

## Elk ------------------
mod1 <- glmmTMB(PredPropMort ~ Age_Class*Predator_Common + Number_of_Relevant_Predators + (1|random),
                betabinomial(link = "logit"),
                weights=Predation,
                data=Elk1)
summary(mod1)

### Plot --------------------------------
tp <- range(Elk1$Number_of_Relevant_Predators)
nd <- expand.grid(Predator_Common = factor(levels(Elk1$Predator_Common), levels = levels(Elk1$Predator_Common)),
                  Age_Class = factor(levels(Elk1$Age_Class), levels = levels(Elk1$Age_Class)),
                  Number_of_Relevant_Predators = seq(tp[1], tp[2], by = 1)) %>%
  mutate(random = levels(Elk1$random)[1],  # Set random to a default value
         Predation = 0) # Set Predation to 0

pred <- predict(mod1, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

ggplot(data=pred, aes(x=Number_of_Relevant_Predators, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  facet_wrap("Predator_Common") +
  mytheme +
  labs(y="Proportion predation",
       x="Number of predators", col="Age class",
       fill="Age class") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F"))
ggsave("Figures/Elk - Model 1.jpg", width = 6, height = 6, dpi = 700)

## Moose -----------
mod1 <- glmmTMB(PredPropMort ~ Age_Class*Predator_Common + Number_of_Relevant_Predators + (1|random),
                betabinomial(link = "logit"),
                weights=Predation,
                data=Moose1)
summary(mod1)

### Plot -------------------
tp <- range(Moose1$Number_of_Relevant_Predators)
nd <- expand.grid(Predator_Common = factor(levels(Moose1$Predator_Common), levels = levels(Moose1$Predator_Common)),
                  Age_Class = factor(levels(Moose1$Age_Class), levels = levels(Moose1$Age_Class)),
                  Number_of_Relevant_Predators = seq(tp[1], tp[2], by = 1)) %>%
  mutate(random = levels(Moose1$random)[1], 
         Predation = 0) 

pred <- predict(mod1, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

ggplot(data=pred, aes(x = Number_of_Relevant_Predators)) + 
 geom_point(aes(y = fit_trans, color = Age_Class), position = position_dodge(0.25)) +
 geom_errorbar(aes(ymin = low_trans, ymax = up_trans, color = Age_Class), width = 0.2, position = position_dodge(0.25)) +
 facet_wrap("Predator_Common") +
 mytheme +
 labs(y = "Proportion predation",
      x = "Number of predators",
      color = "Age class") +
 scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
 scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
 scale_fill_manual(values = c("#1F78B4", "#7FC97F"))
ggsave("Figures/Moose - Model 1.jpg", width = 6, height = 4, dpi = 700)

## Mule deer ------------
mod1 <- glmmTMB(PredPropMort ~ Age_Class*Predator_Common + Number_of_Relevant_Predators + (1|random),
                betabinomial(link = "logit"),
                weights=Predation,
                data=Muledeer1)
summary(mod1)

### Plot -------------------
tp <- range(Muledeer1$Number_of_Relevant_Predators)
nd <- expand.grid(Predator_Common = factor(levels(Muledeer1$Predator_Common), levels = levels(Muledeer1$Predator_Common)),
                  Age_Class = factor(levels(Muledeer1$Age_Class), levels = levels(Muledeer1$Age_Class)),
                  Number_of_Relevant_Predators = seq(tp[1], tp[2], by = 1)) %>%
  mutate(random = levels(Muledeer1$random)[1], 
         Predation = 0) 

pred <- predict(mod1, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

ggplot(data=pred, aes(x=Number_of_Relevant_Predators, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  facet_wrap("Predator_Common") +
  mytheme +
  labs(y="Proportion predation",
       x="Number of predators", col="Age class",
       fill="Age class") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F"))
ggsave("Figures/Muledeer - Model 1.jpg", width = 6, height = 6, dpi = 700)

## Pronghorn ------------
mod1 <- glmmTMB(PredPropMort ~ Age_Class*Predator_Common + Number_of_Relevant_Predators + (1|random),
                betabinomial(link = "logit"),
                weights=Predation,
                data=Pronghorn1)
summary(mod1)
# Analysis was concluded:
#all models have the same number of relevant predators

## Whitetail ----------------
mod1 <- glmmTMB(PredPropMort ~ Age_Class*Predator_Common + Number_of_Relevant_Predators + (1|random),
                betabinomial(link = "logit"),
                weights=Predation,
                data=Whitetail1)
summary(mod1)

### Plot -------------------
tp <- range(Whitetail1$Number_of_Relevant_Predators)
nd <- expand.grid(Predator_Common = factor(levels(Whitetail1$Predator_Common), levels = levels(Whitetail1$Predator_Common)),
                  Age_Class = factor(levels(Whitetail1$Age_Class), levels = levels(Whitetail1$Age_Class)),
                  Number_of_Relevant_Predators = seq(tp[1], tp[2], by = 1)) %>%
  mutate(random = levels(Whitetail1$random)[1], 
         Predation = 0) 

pred <- predict(mod1, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

ggplot(data=pred, aes(x=Number_of_Relevant_Predators, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  facet_wrap("Predator_Common") +
  mytheme +
  labs(y="Predation proportion",
       x="Predator species richness", col="Age class",
       fill="Age class") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F"))
ggsave("Figures/Whitetail - Model 1.jpg", width = 6, height = 6, dpi = 700)

# Models 2-4 ---------------
# Prep
prediction3data <- preddata %>%
  mutate(Percent_Forest_Cover = Percent_Forest_Cover * 100)

# Isolate Ungulates
common_names <- unique(prediction3data$Common_Name)
df_names <- c("Caribou", "Elk", "Moose", "Muledeer", "Whitetail", "Bighorn", "Pronghorn")

split_data <- lapply(common_names, function(name) {
  prediction3data %>%
    filter(Common_Name == name) %>%
    na.omit() %>%
    mutate_at(vars(Age_Class, Predator_Common), as.factor)
})

list2env(setNames(split_data, paste0(df_names, "3")), envir = .GlobalEnv)

## Bighorn --------------
# Analysis was not run because the number of relevant predators was equal at all sites

## Caribou --------------
mod3HFI <- glmmTMB(PredPropMort ~ Human_Footprint_Index * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Caribou3)
summary(mod3HFI)

mod3PFC <- glmmTMB(PredPropMort ~ Percent_Forest_Cover * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Caribou3)
summary(mod3PFC)
mod3TRI <- glmmTMB(PredPropMort ~ Average_Terrain_Ruggedness * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Caribou3)
summary(mod3TRI)

### Plot --------------
tt <- range(Caribou3$Number_of_Relevant_Predators)
tp <- range(Caribou3$Average_Terrain_Ruggedness)

nd <- expand.grid(
  Predator_Common = factor(levels(Caribou3$Predator_Common), levels = levels(Caribou3$Predator_Common)),
  Age_Class = factor(levels(Caribou3$Age_Class), levels = levels(Caribou3$Age_Class)),
  Average_Terrain_Ruggedness = seq(tp[1], tp[2], by = 1),
  Number_of_Relevant_Predators = seq(tt[1], tt[2], by = 1)
) %>%
  mutate(random = levels(Caribou3$random)[1], 
         Predation = 0)  

pred <- predict(mod3TRI, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

ggplot(data=pred, aes(x=Average_Terrain_Ruggedness, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  facet_wrap("Predator_Common") +
  mytheme +
  labs(y="Proportion predation",
       x="TRI", col="Age class",
       fill="Age class") +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F"))
ggsave("Figures/Caribou - Model 3.jpg", width = 6, height = 4, dpi = 700)

# Alternative model without number of predators to make graph nicer
mod3TRI <- glmmTMB(PredPropMort ~ Average_Terrain_Ruggedness * Predator_Common + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Caribou3)

# Plot Prep
tt <- range(Caribou3$Number_of_Relevant_Predators)
tp <- range(Caribou3$Average_Terrain_Ruggedness)

nd <- expand.grid(Predator_Common = factor(levels(Caribou3$Predator_Common), levels = levels(Caribou3$Predator_Common)),
            Age_Class = factor(levels(Caribou3$Age_Class), levels = levels(Caribou3$Age_Class)),
            Average_Terrain_Ruggedness = seq(tp[1], tp[2], by=1),
            Number_of_Relevant_Predators = seq(tt[1], tt[2], by=1)) %>%
  mutate(Predation = 0) 

pred <- predict(mod3TRI, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

ggplot(data=pred, aes(x=Average_Terrain_Ruggedness, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  facet_wrap("Predator_Common") +
  mytheme +
  labs(y="Predation proportion",
       x="TRI", col="Age class",
       fill="Age class") +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F"))
ggsave("Figures/Caribou - Model 3.jpg", width = 6, height = 4, dpi = 700)

## Elk -----------------
mod3HFI <- glmmTMB(PredPropMort ~ Human_Footprint_Index * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Elk3)
summary(mod3HFI)
mod3PFC <- glmmTMB(PredPropMort ~ Percent_Forest_Cover * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Elk3)
summary(mod3PFC)
mod3TRI <- glmmTMB(PredPropMort ~ Average_Terrain_Ruggedness * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Elk3)
summary(mod3TRI)
# Analysis was concluded:
#Model 1 AIC was better than Models 2-4

## Moose ---------------
mod3HFI <- glmmTMB(PredPropMort ~ Human_Footprint_Index * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Moose3)
summary(mod3HFI)
mod3PFC <- glmmTMB(PredPropMort ~ Percent_Forest_Cover * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Moose3)
summary(mod3PFC)
mod3TRI <- glmmTMB(PredPropMort ~ Average_Terrain_Ruggedness * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Moose3)
summary(mod3TRI)
# Analysis was concluded:
#Model 1 AIC was better than Models 2-4

## Mule deer ---------------
mod3HFI <- glmmTMB(PredPropMort ~ Human_Footprint_Index * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Muledeer3)
summary(mod3HFI)
mod3PFC <- glmmTMB(PredPropMort ~ Percent_Forest_Cover * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Muledeer3)
summary(mod3PFC)
mod3TRI <- glmmTMB(PredPropMort ~ Average_Terrain_Ruggedness * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Muledeer3)
summary(mod3TRI)
# Analysis was concluded:
#Model 1 AIC was better than Models 2-4

## Pronghorn ---------------
# Analysis was not run because the number of relevant predators was equal at all sites

## Whitetail -------------------------
mod3HFI <- glmmTMB(PredPropMort ~ Human_Footprint_Index * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Whitetail3)
summary(mod3HFI)
mod3PFC <- glmmTMB(PredPropMort ~ Percent_Forest_Cover * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Whitetail3)
summary(mod3PFC)
mod3TRI <- glmmTMB(PredPropMort ~ Average_Terrain_Ruggedness * Predator_Common + Number_of_Relevant_Predators + Age_Class + (1|random),
                   betabinomial(link = "logit"),
                   weights=Predation,
                   data=Whitetail3)
summary(mod3TRI)
# Analysis was concluded:
#Model 1 AIC was better than Models 2-4

# Model 5 --------------------------------------------------------
prediction2data <- alldata %>% 
  filter(Age_Class != "Mixed", Number_of_Relevant_Predators != 0) %>%
  mutate(Pred2 = Number_of_Relevant_Predators * Number_of_Relevant_Predators)

# Isolate Ungulates
common_names <- unique(prediction2data$Common_Name)
df_names <- c("Bighorn", "Caribou", "Elk", "Moose", "Muledeer", "Pronghorn", "Whitetail")

split_data <- lapply(common_names, function(name) {
  prediction2data %>%
    filter(Common_Name == name) %>%
    mutate(PredPropMort = Predation / Known_Mortalities) %>% 
    mutate_at(vars(Age_Class, random), as.factor)
})

list2env(setNames(split_data, paste0(df_names, "5")), envir = .GlobalEnv)

## Bighorn -------------
# Analysis was not run because predator species richness was equal across all sites

## Caribou -------------
mod5 <- glmmTMB(PredPropMort ~ Number_of_Relevant_Predators*Age_Class + (1|random),
                betabinomial(link = "logit"),
                weights=Known_Mortalities,
                data=Caribou5)
summary(mod5)

### Plot ------------------------
tp <- range(Caribou5$Number_of_Relevant_Predators)
nd <- expand.grid(Age_Class=factor(levels(Caribou5$Age_Class), levels=levels(Caribou5$Age_Class)),
                  Number_of_Relevant_Predators = seq(tp[1], tp[2], by=1)) %>% 
  mutate(random = levels(Caribou5$random)[1],
         Known_Mortalities = 0) 

pred <- predict(mod5, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

Caribou1 <- ggplot(data=pred, aes(x=Number_of_Relevant_Predators, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  mytheme +
  labs(y = "Proportion\npredation",
       x="Number of predators", col="Age class",
fill="Age class") +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F")) +
  ggtitle("Caribou") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

## Elk ---------------
mod5 <- glmmTMB(PredPropMort ~ Number_of_Relevant_Predators*Age_Class + (1|random),
                betabinomial(link = "logit"),
                weights=Known_Mortalities,
                data=Elk5)
summary(mod5)

### Plot ------------------------
tp <- range(Elk5$Number_of_Relevant_Predators)
nd <- expand.grid(Age_Class=factor(levels(Elk5$Age_Class), levels=levels(Elk5$Age_Class)),
                  Number_of_Relevant_Predators = seq(tp[1], tp[2], by=1)) %>% 
  mutate(random = levels(Elk5$random)[1],
         Known_Mortalities = 0) 

pred <- predict(mod5, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

Elk1 <- ggplot(data=pred, aes(x=Number_of_Relevant_Predators, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  mytheme +
  labs(y="Proportion\npredation",       
       x="Number of Predators", col="Age class",
fill="Age class") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F")) +
  ggtitle("Elk") +
  theme(legend.position = "right") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

## Moose -------------------------
mod5 <- glmmTMB(PredPropMort ~ Number_of_Relevant_Predators*Age_Class + (1|random),
                betabinomial(link = "logit"),
                weights=Known_Mortalities,
                data=Moose5)
summary(mod5)
# Analysis was concluded:
#moose model failed to converge

## Mule deer -----------------------------------
mod5 <- glmmTMB(PredPropMort ~ Number_of_Relevant_Predators*Age_Class + (1|random),
                betabinomial(link = "logit"),
                weights=Known_Mortalities,
                data=Muledeer5)
summary(mod5)

### Plot ------------------------
tp <- range(Muledeer5$Number_of_Relevant_Predators)
nd <- expand.grid(Age_Class=factor(levels(Muledeer5$Age_Class), levels=levels(Muledeer5$Age_Class)),
                  Number_of_Relevant_Predators = seq(tp[1], tp[2], by=1)) %>% 
  mutate(random = levels(Muledeer5$random)[1],
         Known_Mortalities = 0) 

pred <- predict(mod5, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

Muledeer1 <- ggplot(data=pred, aes(x=Number_of_Relevant_Predators, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  mytheme +
  labs(y="Proportion\npredation",
       x="Number of Predators", col="Age class",
fill="Age class") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F")) +
  ggtitle("Mule deer") +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

## Pronghorn ------------------------------
# Analysis was not run because predator species richness was equal across all sites

## Whitetail ---------------------
mod5 <- glmmTMB(PredPropMort ~ Number_of_Relevant_Predators*Age_Class + (1|random),
                betabinomial(link = "logit"),
                weights=Known_Mortalities,
                data=Whitetail5)
summary(mod5)

### Plot ------------------------
tp <- range(Whitetail5$Number_of_Relevant_Predators)
nd <- expand.grid(Age_Class=factor(levels(Whitetail5$Age_Class), levels=levels(Whitetail5$Age_Class)),
                  Number_of_Relevant_Predators = seq(tp[1], tp[2], by=1)) %>% 
  mutate(random = levels(Whitetail5$random)[1],
         Known_Mortalities = 0) 

pred <- predict(mod5, newdata = nd, se.fit = TRUE, re.form = NA) %>%
  as.data.frame() %>%
  mutate(lower = fit - 1.96 * se.fit,
         upper = fit + 1.96 * se.fit,
         fit_trans = plogis(fit),
         low_trans = plogis(lower),
         up_trans = plogis(upper)) %>%
  bind_cols(nd)

Whitetail1 <- ggplot(data=pred, aes(x=Number_of_Relevant_Predators, y=fit_trans)) +
  geom_ribbon(aes(ymin=low_trans, ymax=up_trans, fill=Age_Class), alpha=0.2) +
  geom_line(aes(col=Age_Class), linewidth = 0.75) + 
  mytheme +
  labs(x = "Number of predators", y="Proportion\npredation", col="Age class",
fill="Age class") +
  scale_x_continuous(breaks=c(0,1,2,3,4,5)) +
  scale_color_manual(values = c("#1F78B4", "#7FC97F")) +
  scale_fill_manual(values = c("#1F78B4", "#7FC97F")) +
  ggtitle("White-tailed deer") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())

# Combination figure
Combo1 <- ggarrange(Caribou1 + rremove("ylab") + rremove("xlab"), 
                    Elk1 + rremove("ylab") + rremove("xlab"), 
                    Muledeer1 + rremove("ylab") + rremove("xlab"),
                    Whitetail1 + rremove("ylab") + rremove("xlab"),
                    labels = NULL,
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "right",
                    align = "hv", 
                    font.label = list(size = 12, color = "black", family = NULL, position = "top"))
Combo1 <- annotate_figure(Combo1, left = textGrob("Predation proportion", rot = 90, vjust = 1, gp = gpar(cex = 1.3)), bottom = textGrob("Predator species richness", gp = gpar(cex = 1.3)))
Combo1
ggsave("Figures/All Species - Model 5.jpg", width = 6.5, height = 6.5, dpi = 700)
