#########################################################################################################################################
# Manuscript title: Perinatal Oxycodone Exposure Causes Long Term Sex-Dependent Changes in Sensory and Reward Processing in Adult Mice  #
# Manuscript authors: Minakova, E. & Mikati, M. et al.                                                                                  #
# Authors: Justin Baldwin and Marwa Mikati                                                                                              #
# 2022-02-10                                                                                                                            #
#########################################################################################################################################

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(viridis)
library(dplyr)
library(lmerTest)
library(bbmle)
lines<-theme(axis.line.x.bottom = element_line(), axis.line.y.left = element_line())

################################ 1. Analysis section: FED3 operant conditioning ################################
#load data

m<-read.csv("LongFormat_FED_Data_NAS.csv", skip=1)
m_long<-pivot_longer(m %>% select(1:22), cols = 8:22, names_to = "day", values_to = "active_nosepoke") %>%
  mutate(day = as.numeric(gsub("Active_Nosepokes_Day", "", day))) %>%
  mutate(ID_day=paste0(ID, "_", day))
m_long<-left_join(m_long, 
                  pivot_longer(m %>% select(c(2, 23:37)), cols = 2:16, names_to = "day", values_to = "inactive_nosepoke")  %>%
                    mutate(day = as.numeric(gsub("Inactive_Nosepokes_Day", "", day))) %>%
                    mutate(ID_day=paste0(ID, "_", day)) %>%
                    select(c(ID_day, inactive_nosepoke)), 
                  by="ID_day")
m_long<-left_join(m_long, 
                  pivot_longer(m %>% select(c(2, 38:52)), cols = 2:16, names_to = "day", values_to = "earned_pellets")  %>%
                    mutate(day = as.numeric(gsub("Pellets_Earned_Day", "", day))) %>%
                    mutate(ID_day=paste0(ID, "_", day)) %>%
                    select(c(ID_day, earned_pellets)), 
                  by="ID_day")
m_long<-left_join(m_long, 
                  pivot_longer(m %>% select(c(2, 53:67)), cols = 2:16, names_to = "day", values_to = "consumed_pellets")  %>%
                    mutate(day = as.numeric(gsub("Pellets_Consumed_Day", "", day))) %>%
                    mutate(ID_day=paste0(ID, "_", day)) %>%
                    select(c(ID_day, consumed_pellets)), 
                  by="ID_day")


################################ 1. Analysis section: FED3 operant conditioning ################################
############# 1.1. Pellets consumed ############################################################################

#does it matter if we use one or the other?
plot(earned_pellets~consumed_pellets, data=m_long)
#shouldnt matter, highly correlated
#cor.test(m_long$earned_pellets, m_long$consumed_pellets)# t = 65.01, df = 679, p-value < 2.2e-16, cor = 0.9282133 

#visual exploration
#Group is the combination of drug (Oxy/Vehicle) and duration (short x long)
m_long$drug<-NA
m_long$duration<-NA

m_long$drug[m_long$Group %in% c("A", "C")]<-"oxy"
m_long$drug[m_long$Group%in% c("B", "D")]<-"vehicle"
m_long$duration[m_long$Group %in% c("A", "B")]<-"short"
m_long$duration[m_long$Group %in% c("C", "D")]<-"long"

m_long$group_unblinded<-paste0(m_long$drug, "_", m_long$duration)

ggplot(m_long, aes(x=day, y=consumed_pellets, group=ID, color=Sex))+
  geom_line(alpha=0.5, size=2, position=position_jitter())+
  theme_tufte()+ylab("Consumed Pellets")+xlab("Day")+
  facet_wrap(~group_unblinded, ncol=4)+lines

ggplot(m_long, aes(x=day, y=consumed_pellets, group=ID, color=Sex))+
  geom_line(alpha=0.3, position=position_jitter())+
  theme_tufte()+ylab("Consumed Pellets")+xlab("Day")+
  facet_wrap(~Sex, ncol=4)+lines

ggplot(m_long %>% mutate(sex_group=paste0(Sex, Group)), aes(x=day, y=consumed_pellets, group=ID, color=sex_group))+
  geom_line(alpha=0.3, position=position_jitter())+
  theme_tufte()+ylab("Consumed Pellets")+xlab("Day")+
  facet_grid(group_unblinded~Sex)+lines

ggplot(m_long %>% mutate(sex_group=paste0(Sex, Group)), aes(x=day, y=consumed_pellets, group=ID, color=sex_group))+
  geom_line(alpha=0.3, position=position_jitter())+
  theme_tufte()+ylab("Consumed Pellets")+xlab("Day")+
  facet_grid(Sex~group_unblinded)+lines


hist(m_long$consumed_pellets) #this looks mostly normal, although with a theoretical lower bound at 0 (mice can't consume negative pellets).
#linear models are probably justified

#backwards model simplification using linear mixed models. ID is random effect for Mouse ID (repeated measurements of individual mice over time). 
#model simplification is performed by eliminating non-significant terms
#model comparisons are done with maximum likelihood, the final model is refitted with restricted maximum likelihood


#lmm 
m1_p<-lmer(consumed_pellets~day + Group + Sex +day:Group + day:Sex + Group:Sex + 
             day:Group:Sex + (1|ID), 
           data = m_long, REML = F)
m2_p<-lmer(consumed_pellets~day + Group + Sex +day:Group + day:Sex + Group:Sex + 
              (1|ID), 
            data = m_long, REML = F)
anova(m1_p, m2_p, method="LRT")#this is default method
# Data: m_long
# Models:
#   m2_pm: consumed_pellets ~ day + Group + Sex + day:Group + day:Sex + Group:Sex + (1 | ID)
# m1_p: consumed_pellets ~ day + Group + Sex + day:Group + day:Sex + Group:Sex + day:Group:Sex + (1 | ID)
# npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# m2_pm   15 3139.1 3206.9 -1554.5   3109.1                    
# m1_p    18 3139.0 3220.4 -1551.5   3103.0 6.098  3     0.1069

m3_p<-lmer(consumed_pellets~day + Group + Sex +day:Group + day:Sex + 
             (1|ID), 
           data = m_long, REML = F)
anova(m2_p, m3_p) #no group x sex
# > anova(m2_p, m3_p) #no group x sex
# Data: m_long
# Models:
#   m3_p: consumed_pellets ~ day + Group + Sex + day:Group + day:Sex + (1 | ID)
# m2_p: consumed_pellets ~ day + Group + Sex + day:Group + day:Sex + Group:Sex + (1 | ID)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# m3_p   12 3136.2 3190.5 -1556.1   3112.2                     
# m2_p   15 3139.1 3206.9 -1554.5   3109.1 3.1284  3     0.3722

m4_p<-lmer(consumed_pellets~day + Group + Sex +day:Group + 
             (1|ID), 
           data = m_long, REML = F)
anova(m3_p, m4_p) #marginal day x sex
# Data: m_long
# Models:
#   m4_p: consumed_pellets ~ day + Group + Sex + day:Group + (1 | ID)
# m3_p: consumed_pellets ~ day + Group + Sex + day:Group + day:Sex + (1 | ID)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# m4_p   11 3137.4 3187.1 -1557.7   3115.4                       
# m3_p   12 3136.2 3190.5 -1556.1   3112.2 3.1648  1    0.07524 .
# -


m5_p<-lmer(consumed_pellets~day + Group + Sex +
             (1|ID), 
           data = m_long, REML = F)
anova(m4_p, m5_p)
# Data: m_long
# Models:
#   m5_p: consumed_pellets ~ day + Group + Sex + (1 | ID)
# m4_p: consumed_pellets ~ day + Group + Sex + day:Group + (1 | ID)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# m5_p    8 3133.6 3169.7 -1558.8   3117.6                     
# m4_p   11 3137.4 3187.1 -1557.7   3115.4 2.1678  3     0.5383

m6_p<-lmer(consumed_pellets~day + Group  +(1|ID), 
           data = m_long, REML = F)
anova(m5_p, m6_p) # there is a sig effect of sex
# Data: m_long
# Models:
#   m6_p: consumed_pellets ~ day + Group + (1 | ID)
# m5_p: consumed_pellets ~ day + Group + Sex + (1 | ID)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# m6_p    7 3142.2 3173.9 -1564.1   3128.2                        
# m5_p    8 3133.6 3169.7 -1558.8   3117.6 10.669  1   0.001089 **

m7_p<-lmer(consumed_pellets~day + Sex  +
             (1|ID), 
           data = m_long, REML = F)
anova(m5_p, m7_p) # there is there is no effect of group
# Data: m_long
# Models:
#   m7_p: consumed_pellets ~ day + Sex + (1 | ID)
# m5_p: consumed_pellets ~ day + Group + Sex + (1 | ID)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# m7_p    5 3127.9 3150.5 -1558.9   3117.9                     
# m5_p    8 3133.6 3169.7 -1558.8   3117.6 0.3178  3     0.9567
m8_p<-lmer(consumed_pellets~ Sex  +
              
              (1|ID), 
            data = m_long, REML = F)
anova(m7_p, m8_p) # there is there is an effect of day

#final model includes day and sex


#just be certain the marginal interaction between sex and day didnt get prematurely discarded, we can add it back in and give it one last shot

m9_p<-lmer(consumed_pellets~day + Sex +Sex:day + (1|ID), 
           data = m_long, REML = F)
anova(m7_p, m9_p) # there is there is an effect of day
#p0.0761 Nope

#refit final model with reml

m7_p_final<-lmer(consumed_pellets~day + Sex  +(1|ID), 
                 data = m_long, REML = T)
summary(m7_p_final)
plot(m7_p_final)#one large outlier, and perhaps a mild fan shape from low to central fitted values,but then remains smudgy.
#this might be a case of very mild heteroscedasticity. proceed with caution.

############# 1.2. Nosepoke accuracy ############################################################################

m_long$acc_nosepoke<-m_long$active_nosepoke/c(m_long$active_nosepoke+m_long$inactive_nosepoke)

id_final_stats_nosepoke_acc<-m_long %>% 
  group_by(ID, Sex, Group) %>%  
  filter(!is.na(acc_nosepoke)) %>%
  summarise(final_acc=acc_nosepoke[day==Total.days.of.training], 
            day=max(day))
id_final_stats_nosepoke_acc
#Mean of accuracy of nosepoke and mean plot 

MeanAccuracy <- m_long %>% group_by(m_long$day,m_long$group_unblinded) %>% summarize(mean_accuracy=mean(acc_nosepoke, na.rm=TRUE))
MeanAccuracy <- as.tibble(MeanAccuracy)
MeanAccuracy <- rename(MeanAccuracy, 'Day' = 'm_long$day', 'Group' = 'm_long$group_unblinded')

ggplot(MeanAccuracy, aes(x=Day, y=mean_accuracy, group=Group, color=Group))+
         geom_line(alpha=0.3, position=position_jitter())+
         theme_tufte()+ylab("Mean Accuracy")+xlab("Day")+
         facet_wrap(~Group, ncol=1)+lines


ggplot(MeanAccuracy, aes(x=Day, y=mean_accuracy, group=Group, color=Group))+
         geom_line(alpha=0.3, position=position_jitter())+
         theme_tufte()+ylab("Mean Accuracy")+xlab("Day")+
         facet_wrap(~Group)+lines

#sex and accuracy # adding final accuracy at each mouses end of trial
ggplot(m_long, aes(x=day, y=acc_nosepoke, group=ID, color=Sex))+
         geom_line(alpha=0.3, position = position_jitter())+
         geom_point(data = id_final_stats_nosepoke_acc,
                    aes(x=day, y=final_acc, color=Sex), size=2)+
         scale_color_manual(values=c("M"='red', "F"='blue'))+theme_tufte()+
         ylab("Accuracy Nosepoke")+xlab("Day")+
         facet_wrap(~Sex)+lines


#group and accuracy
ggplot(m_long, aes(x=day, y=acc_nosepoke, group=ID, color=Group))+
         geom_line(alpha=0.5, position = position_jitter())+
         geom_point(data = id_final_stats_nosepoke_acc,
                    aes(x=day, y=final_acc, color=Group), size=2, alpha=0.8 ,position = position_jitter())+
         theme_tufte()+
         ylab("Accuracy Nosepoke")+xlab("Day")+
         facet_wrap(~Group, ncol = 4)

ggplot(m_long, aes(x=day, y=acc_nosepoke, group=ID, color=Group))+
         geom_line(alpha=0.5, position = position_jitter())+
         geom_point(data = id_final_stats_nosepoke_acc,
                    aes(x=day, y=final_acc, color=Group), size=2, alpha=0.8 ,position = position_jitter())+
         theme_tufte()+
         ylab("Accuracy Nosepoke")+xlab("Day")+
         facet_wrap(~Group, ncol = 1)

#sex x group and accuracy
ggplot(m_long, aes(x=day, y=acc_nosepoke, group=ID, color=Sex))+
         geom_line(alpha=0.5, position = position_jitter())+
         geom_point(data = id_final_stats_nosepoke_acc,
                    aes(x=day, y=final_acc, color=Sex), size=2, alpha=0.8 ,position = position_jitter())+
         theme_tufte()+
         ylab("Accuracy Nosepoke")+xlab("Day")+
         facet_wrap(~Group)

m %>% count(Sex, Group)
m %>% count(Sex)

hist(m_long$acc_nosepoke)
#this looks relatively normal, event though accuracy is bounded by 1 and 0.

#an lmm might be tolerable in this case, unless we want to consider other options.......

m1_n<-lmer(acc_nosepoke~day + Group + Sex +day:Group + day:Sex + Group:Sex + day:Group:Sex + (1|ID), 
           data = m_long, REML = F)
m2_n<-lmer(acc_nosepoke~day + Group + Sex +day:Group + day:Sex + Group:Sex+ (1|ID), 
           data = m_long, REML = F)
anova(m1_n, m2_n) #no three way
# Data: m_long
# Models:
#   m2_n: acc_nosepoke ~ day + Group + Sex + day:Group + day:Sex + Group:Sex + (1 | ID)
# m1_n: acc_nosepoke ~ day + Group + Sex + day:Group + day:Sex + Group:Sex + day:Group:Sex + (1 | ID)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
# m2_n   15 -828.17 -760.32 429.09  -858.17                     
# m1_n   18 -825.12 -743.69 430.56  -861.12 2.9462  3        0.4

m3_n<-lmer(acc_nosepoke~day + Group + Sex +day:Group +day:Sex+  (1|ID), 
           data = m_long, REML = F)
anova(m2_n, m3_n)  #marginal group x sex
# Data: m_long
# Models:
#   m3_n: acc_nosepoke ~ day + Group + Sex + day:Group + day:Sex + (1 | ID)
# m2_n: acc_nosepoke ~ day + Group + Sex + day:Group + day:Sex + Group:Sex + (1 | ID)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)  
# m3_n   12 -826.62 -772.34 425.31  -850.62                       
# m2_n   15 -828.17 -760.32 429.09  -858.17 7.5513  3    0.05625 .

m4_n<-lmer(acc_nosepoke~day + Group + Sex +day:Group +(1|ID), data = m_long, REML = F)
anova(m3_n, m4_n)  #no day x sex
# Data: m_long
# Models:
#   m4_n: acc_nosepoke ~ day + Group + Sex + day:Group + (1 | ID)
# m3_n: acc_nosepoke ~ day + Group + Sex + day:Group + day:Sex + (1 | ID)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
# m4_n   11 -828.13 -778.37 425.07  -850.13                     
# m3_n   12 -826.62 -772.34 425.31  -850.62 0.4883  1     0.4847

m5_n<-lmer(acc_nosepoke~day + Group + Sex + (1|ID), data = m_long, REML = F)
anova(m4_n, m5_n)  #no day x group
# Data: m_long
# Models:
#   m5_n: acc_nosepoke ~ day + Group + Sex + (1 | ID)
# m4_n: acc_nosepoke ~ day + Group + Sex + day:Group + (1 | ID)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
# m5_n    8 -831.06 -794.87 423.53  -847.06                     
# m4_n   11 -828.13 -778.37 425.07  -850.13 3.0691  3     0.3811


m6_n<-lmer(acc_nosepoke~day + Group  + (1|ID), data = m_long, REML = F)
anova(m5_n, m6_n)  #marginal sex effect

#just to be sure we can add that previously marginal interaction between group and sex back in
m5_n_2<-lmer(acc_nosepoke~day + Group + Sex+Group:Sex + (1|ID), data = m_long, REML = F)
anova(m5_n, m5_n_2)  #it's still marginal, no dice. proceed simplifying, but be cautious
# ata: m_long
# Models:
#   m5_n: acc_nosepoke ~ day + Group + Sex + (1 | ID)
# m5_n_2: acc_nosepoke ~ day + Group + Sex + Group:Sex + (1 | ID)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)  
# m5_n      8 -831.06 -794.87 423.53  -847.06                       
# m5_n_2   11 -832.73 -782.97 427.36  -854.73 7.6648  3    0.05347 .

m7_n<-lmer(acc_nosepoke~day   + (1|ID), data = m_long, REML = F)
anova(m7_n, m6_n)  #no effect of group
# Data: m_long
# Models:
#   m7_n: acc_nosepoke ~ day + (1 | ID)
# m6_n: acc_nosepoke ~ day + Group + (1 | ID)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
# m7_n    4 -834.49 -816.39 421.24  -842.49                     
# m6_n    7 -830.10 -798.43 422.05  -844.10 1.6109  3     0.6569



m8_n<-lmer(acc_nosepoke~1   + (1|ID), data = m_long, REML = F)
anova(m7_n, m8_n)  #day has an effect
# Data: m_long
# Models:
#   m8_n: acc_nosepoke ~ 1 + (1 | ID)
# m7_n: acc_nosepoke ~ day + (1 | ID)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)  
# m8_n    3 -832.09 -818.51 419.04  -838.09                       
# m7_n    4 -834.49 -816.39 421.24  -842.49 4.3996  1    0.03595 *


m7_n_final<-lmer(acc_nosepoke~day   + (1|ID), data = m_long, REML = T)
summary(m7_n_final)
plot(m7_n_final)#no obvious heteroscedasticity

############# 1.3. Number of animals meeting criterion ##########################################################

m$Met.Criteria_num<-ifelse(m$Met.Criteria=="Yes", 1, 0)
m$Days.to.criterion_num<-as.numeric(gsub("Days.to.criterion", "", m$Days.to.criterion))

m %>% filter(!is.na(Days.to.criterion_num)) %>%  count(Sex, Group)
m %>%  count(Sex, Group)

mq1<-glm(Met.Criteria_num~Sex*Group , 
         data=m , 
         family = "binomial")
mq2<-glm(Met.Criteria_num~Sex+Group , 
         data=m , 
         family = "binomial")
mq3<-glm(Met.Criteria_num~Sex , 
         data=m , 
         family = "binomial")
mq4<-glm(Met.Criteria_num~Group , 
         data=m , 
         family = "binomial")
mq5<-glm(Met.Criteria_num~1 , 
         data=m , 
         family = "binomial")

#test of interaction
anova(mq1, mq2, test="LRT")
# Analysis of Deviance Table
# 
# Model 1: Met.Criteria_num ~ Sex * Group
# Model 2: Met.Criteria_num ~ Sex + Group
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1        49     59.837                     
# 2        52     62.112 -3  -2.2746   0.5174

#test of group
anova(mq2, mq3, test="LRT")
# Analysis of Deviance Table
# 
# Model 1: Met.Criteria_num ~ Sex + Group
# Model 2: Met.Criteria_num ~ Sex
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1        52     62.112                     
# 2        55     65.401 -3  -3.2887   0.3492



#nope
#test of sex
anova(mq5, mq3, test="LRT")
# Analysis of Deviance Table
# 
# Model 1: Met.Criteria_num ~ 1
# Model 2: Met.Criteria_num ~ Sex
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1        56     69.468                       
# 2        55     65.401  1   4.0677  0.04371 *

summary(mq3)
plot(mq3)


ggplot(m %>% group_by(Sex) %>% count(Met.Criteria_num) %>% mutate(p=n/sum(n)), 
              aes(x=Sex, fill=factor(Met.Criteria_num, levels=c("1", "0")), y=p))+
         geom_col(position = position_stack())+
         scale_fill_manual(name="Sex", 
                           values = c("0"="red", "1"="blue"), 
                           labels=rev(c("Fail", "Succeed")))+ylab("Porportion")+theme_tufte()

       
ggplot(m %>% group_by(Group) %>% count(Met.Criteria_num) %>% mutate(p=n/sum(n)), 
              aes(x=Group, fill=factor(Met.Criteria_num, levels=c("1", "0")), y=p))+
         geom_col(position = position_stack())+
         scale_fill_manual(name="Sex", 
                           values = c("0"="red", "1"="blue"), 
                           labels=rev(c("Fail", "Succeed")))+ylab("Porportion")+theme_tufte()



##########################################################################################
################################ 2. Analysis section: CPP ################################
##########################################################################################


K_long<-read.csv("CPP.csv")
names(K_long)[1]<-"Group"
#str(K_long)
K_long$Group<-as.character(K_long$Group) #must coerce to correct object class (factor/character NOT interger/numeric), this keeps the variable a discrete variable

ggplot(K_long, aes(x=factor(Group), y=Preference,color=Sex))+
  geom_violin(aes(x=factor(Group), y=Preference,color=Sex), position = position_dodge(width=0.5))+
  geom_point(aes(x=factor(Group), y=Preference,color=Sex), position = position_dodge(width=0.5),alpha=0.3)+
  scale_color_manual(values=c("M"='red', "F"='blue'))+theme_tufte()+
  ylab("Preference")+xlab("Group")+lines


hist(K_long$Preference)
#no bounds, normal enough





#lm on drug, duration, and sex 

m1_p<-lm(Preference~  Drug + Sex + Duration + Drug:Sex + Duration:Sex + Drug:Duration + Drug:Duration:Sex, data = K_long)
m2_p<-lm(Preference~  Drug + Sex + Duration + Drug:Sex + Duration:Sex + Drug:Duration, data = K_long)
anova(m1_p,m2_p, test="F")
#No significant 3-way interaction, F=2.2387 p=0.1444
summary(m1_p)

m3_p<-lm(Preference~ Drug + Sex + Duration + Drug:Sex + Duration:Sex, data = K_long)
anova(m2_p,m3_p, test="F")
summary(m2_p)

#Drug:Duration is marginally significant F=3.1214 p=0.08652

m4_p<-lm(Preference~ Drug + Sex + Duration + Drug:Sex + Drug:Duration, data = K_long)
anova(m4_p,m2_p, test="F") 

#no effect of duration:Sex F=0.1734 p=0.6798

m5_p<-lm(Preference~ Drug + Sex + Duration + Duration:Sex + Drug:Duration, data = K_long)

anova(m5_p,m2_p) 

#no effect of Drug:Sex F=2.7119 p=0.1091

#We had the marginal drug:duration

 m6_p<-lm(Preference~ Drug + Sex + Duration+ Drug:Duration, data = K_long)
 m7_p<-lm(Preference~ Drug + Sex + Duration, data = K_long)
# anova(m6_p,m7_p) still marginal
summary(m6_p)


m8_p<-lm(Preference~ Drug + Duration, data = K_long)
anova(m7_p,m8_p)
#significant effect of sex F=10.309  p=0.002786 
summary(m8_p)


m9_p<-lm(Preference~ Sex+Drug, data = K_long)
anova(m9_p,m7_p)
#no effect of duration, F=1.9621  , P = 0.1698




m10_p<-lm(Preference~ Sex+Duration, data = K_long)
anova(m10_p,m7_p)
#no effect of drug, F=0.2357   , P = 0.6302

#######################
########End###########


