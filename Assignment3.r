IrrigationExperimentlong <- read.csv("IrrigationExperimentlong.csv",header = FALSE)
IrrigationExperimentlong[1:3] -> IrrigationExperimentlong
#NB: Data already converted to long format in excel
colnames(IrrigationExperimentlong) <- c("Fertiliser","Block","Yield")
IrrigationExperimentlong$Block<-as.factor(IrrigationExperimentlong$Block)

anova(aov(data = IrrigationExperimentlong,Yield~Block + Fertiliser))
plot(aov(data = IrrigationExperimentlong,Yield~Fertiliser+Block))

#Calculating Relative Efficiency - Irrigation

b <-  8 #blocks
a <-  6 #factors

ms_blocks <- 62930

sigma_r <- (((b-1)*ms_blocks)+(b*(a-1)*sigma_b))/((a*b)-1)
sigma_b <-  4546 #equal to ms_error

dfr <- 42
dfb <- 35

R_irrigation <- (((dfb + 1)*(dfr + 3)*sigma_r)/((dfb + 3)*(dfr + 1)*(sigma_b)))
R_irrigation

#Comparing two treatment means
n <- 8
u_trickle <- mean(IrrigationExperimentlong$Yield[IrrigationExperimentlong$Fertiliser == "Trickle"])
u_basin <- mean(IrrigationExperimentlong$Yield[IrrigationExperimentlong$Fertiliser == "Basin"])
sd_trickle <- sd(IrrigationExperimentlong$Yield[IrrigationExperimentlong$Fertiliser == "Trickle"])
sd_basin <- sd(IrrigationExperimentlong$Yield[IrrigationExperimentlong$Fertiliser == "Basin"])
se_ut <- sqrt((sd_trickle^2/n)+(sd_basin^2/n))
se_ut

redlight <- read.csv("redlightexperiment.csv")
redlight$Intersection <- as.factor(redlight$Intersection)
redlight$Time.Period <- as.factor(redlight$Time.Period)

# ANOVAs for relative efficency R for time period and intersection

timeperiod_anova <- anova(aov(data = redlight,Unused.Red.Light.Time~Red.Light.Sequence + Time.Period))
intersection_anova <- anova(aov(data = redlight,Unused.Red.Light.Time~Intersection + Red.Light.Sequence))

#Calculating Relative Efficiency - time period

b <-  5 #blocks
a <-  5 #factors

ms_blocks <- 272.917 #from ANOVA

sigma_b <-  5.520 #equal to ms_error
sigma_r <- (((b-1)*ms_blocks)+(b*(a-1)*sigma_b))/((a*b)-1)


dfr <- 20 #from CRD anova
dfb <- 16 #from two-way anova with blocking

R_timeperiod <- (((dfb + 1)*(dfr + 3)*sigma_r)/((dfb + 3)*(dfr + 1)*(sigma_b)))
R_timeperiod

#Calculating Relative Efficiency - intersection

b <-  5 #blocks
a <-  5 #factors

ms_blocks <- 4.557 #from ANOVA

sigma_b <-  72.610  #equal to ms_error
sigma_r <- (((b-1)*ms_blocks)+(b*(a-1)*sigma_b))/((a*b)-1)


dfr <- 20 #from CRD anova
dfb <- 16 #from two-way anova with blocking

R_intersection <- (((dfb + 1)*(dfr + 3)*sigma_r)/((dfb + 3)*(dfr + 1)*(sigma_b)))
R_intersection