
### some useful libraries    ######


        library(nlme)     # includes lme()
        library(lme4)    # includes lmer(), glmer()
        library(car)    # includes Anova()
        library(coin)   # includes Fisher Pitman Test           
        library(gplots) # necessary for plotmeans()
        library(LMERConvenienceFunctions)  ## provides different p-value calculations for lmer()
        library(pgirmess) # includes PermTest()

        library(sjstats)   # overdisp()
        
######### Useful code for model diagnostics of mixed-effects models ############

model=a1 ## enter the name of the model you want to diagnose

par(mfrow=c(1,3)) 
# TESTs FOR NORMAL DISTRIBUTION for lmer  #
qqPlot(residuals(model), main="QQ Plot", ylab="Residuals")
hist(residuals(model), main="histogram", ylab="Residuals")

##### HOMOGENEITY OF VARIANCE for lmer -----------------------#######
fitted.values=fitted(model)
residual.values=residuals(model)
plot(residual.values~fitted.values, main="Residuals vs. Fitted values")
abline(mean(residual.values),0)
##### -----------------------####### 


####################### Statistics - Tlaxcala 21.11.2019 ##################

## example kittens.short - LM.LMM.GLMM_M1.M2_2019__B.xlsx

citation("coin")  #para imprimir la cita de un paquete.

data1=read.table(file.choose(),header=TRUE,sep=",", dec=".") # read the data into R

str(data1)

data1$x1

#oneway_test -> Fischer Pitman test Equivalent to T-test for in and dependent variables
# (agonistic.behaviour~sex1, -> dependent variable ~ independent variable,
# distribution=approximate(nresample=10000), -> tipo de permutacion, approximate es montecarlo, exact lo calcula de tu sample size
# data1 -> tu archivo de datos

oneway_test (agonistic.behaviour~sex1,
             distribution=approximate(nresample=10000), data1)

plotmeans(agonistic.behaviour~sex1, data1, connect=F) # plot la medias con 50% confidence

#conclusion: The agonistic behaviour of females is significantly higher than in males
#            Fischer-Pitman test: F=, P<.

# factor(time) and factor(subjetcs) para tranformar un Int (numero entero) a Factor
# |factor(subjects), esta es una forma the considerar la identidad de los individuos (parear) 
# el análisis -> ESTRATIFICAR

oneway_test (bodymass~factor(time)|factor(subjects),
             distribution=approximate(nresample=10000), data1)

plotmeans(bodymass~time, data1, connect=F)

######## test for correlation between continuous variables #######
###########linear regression model parametric ##########
plot(y2~x2, data1)

a1=lm(y2~x2,data1)

# visualizar la distribución y otras características de los datos del modelo
plot(a1)
summary(a1)    # para extraer el valor de F y R2, pero no  para obtener los p values
anova(a1)      # anova con a minúscula NO RECOMENDABLE especialmente con datos no balanceados
Anova(a1, type="2")   # para calcular el valor de p
                      # type=2 cuando datos no balanceados y type=3 cuando tenemos interacciones

#conclusion: there was a significant and negative correlation between the frequency of agonstic behavouir and
# the cortisol concentration (linear regression: F1,26= 129.1, R2 = 0.832, P < 0.001) 

######### Useful code for model diagnostics of mixed-effects models ############

model=a1 ## enter the name of the model you want to diagnose

par(mfrow=c(1,3)) 
# TESTs FOR NORMAL DISTRIBUTION for lmer  #
qqPlot(residuals(model), main="QQ Plot", ylab="Residuals")
hist(residuals(model), main="histogram", ylab="Residuals")

#___________________________________________

##### HOMOGENEITY OF VARIANCE for lmer -----------------------#######
fitted.values=fitted(model)
residual.values=residuals(model)
plot(residual.values~fitted.values, main="Residuals vs. Fitted values")
abline(mean(residual.values),0)
##### -----------------------#######  


### Permutational test for correlation analysis

a1=lm(y2~x2,data1) # necesitamos correr el modelo linear también
summary(a1)        # para extraer el valor de R2
PermTest(a1, B=1000)   # para calcular la p value con permutations

#conclusion: there was a significant and negative correlation between the frequency of agonstic behavouir and
# the cortisol concentration (linear regression with permutation: R2 = 0.832, P < 0.001)

############Comparison between means ###################
## Three group comparison

a2=lm(behaviour.b~group, data1) # corremos el modelo linear 
summary(a2)
PermTest(a2, B=1000) # obtenemos la p value con permutations 

plotmeans(behaviour.b~group, data1, connect=F)

#conclusion: there was a significant difference between the three groups (linear regression with permutation:P < 0.05)

### Now the pos hoc analysis
# group1 vs group 2
a2.1=lm(behaviour.b~group, 
      subset(data1, group=="group1"| group=="group2")) # corremos el modelo linear 
PermTest(a2.1, B=1000)

# group2 vs group 3
a2.2=lm(behaviour.b~group, 
        subset(data1, group=="group2"| group=="group3")) # corremos el modelo linear
PermTest(a2.2, B=1000)

# group1 vs group 3
a2.3=lm(behaviour.b~group, 
        subset(data1, group=="group1"| group=="group3")) # corremos el modelo linear
PermTest(a2.3, B=1000)

############# Linear mixed effects ####
# Acepta missing values
|   or
== equals to
!= unquals to

plotmeans(behaviour.a~timestep, data1)  # 

b1=lme(behaviour.a~timestep, random=~1 |ind,   # random es para considerar que hay mediciones repetidas, i.e. corregir dependencias. Es como el anova de medias repetidas
       subset (data1, behaviour.a != "NA"))    # esto es para indicar que no considere los NA
                                                # porque lme si es sensible a missing data
PermTest(b1, B=1000)

### Now the pos hoc analysis
# group1 vs group 2
b1.1=lme(behaviour.a~timestep, random=~1 |ind, 
        subset(data1, behaviour.a != "NA"
               & (timestep=="time1"| timestep=="time2"))) # corremos el modelo linear 

PermTest(b1.1, B=1000)

# group2 vs group 3
b1.2=lme(behaviour.a~timestep, random=~1 |ind,
        subset(data1, (timestep=="time2"| timestep=="time3"))) # corremos el modelo linear
PermTest(b1.2, B=1000)

# group1 vs group 3
b1.3=lme(behaviour.a~timestep, random=~1 |ind,
        subset(data1, behaviour.a != "NA"
               & (timestep=="time1"| timestep=="time3"))) # corremos el modelo linear 

PermTest(b1.3, B=1000)




# en el escenario hipotético que tenemos este experimento:
#n= 30 individuals
#behaviour = dependent variable 
#time = 3 timesteps / factor()
#ind= individual identity (random factor)
#group = always 6 individuals in a cage

#library(nlme)
#library(pgirmess)
# a1=lme(behaviour~time,
# random=~1|group/ind)
# PermTest(a1, B= 1000)
#
lmer()  # recibe todos los random factor que sean, pero no acepta permutaciones,
# entonces las p values se deben sacar con Anova
#library (lme2) 
#a1=lmer(behaviour~factor(time) +
# (1|ind) + 1|group
# (1|group), data1)
# Anova(a1)


############ Multifactorial statistics

data2=read.table(file.choose(),header=TRUE,sep=",", dec=".") # read the data into R

# to test colinearity between predictors to diagnose this prerrquisite
# variance inflation factors < 5 is ok
a1=lme(longevity~
               thorax+
               density+
               #thorax*density+     # * indica interacción
               activity,
       random=~1|cage, data2)

vif.mer(a1)    # obtener la función de vif.for.glme.r archive

a1=lme(log(longevity+1)~
               thorax+
               density+
               thorax*density+     # * indica interacción
               activity,
       random=~1|cage, data2)

# aquí podemos correr las gráficas para diagnosticar si cumplen los prerrequisitos de normalidad
# como no los cumple, podemos transformar log longevity. Ahora podemos hacer un paramétrico,
# o hacer el permutational test

PermTest(a1, B=1000)

par(mfrow=c(1,2))
plotmeans(longevity~density, main= "long thorax", subset(data2, thorax=="long" ))
plotmeans(longevity~density,  main= "short thorax",subset(data2, thorax=="short" ))


# for long thorax

a1.1=lme(log(longevity+1)~
               density,
               random=~1|cage, 
       subset(data2, thorax=="long" ))

PermTest(a1.1, B=1000)

# for short thorax

a1.2=lme(log(longevity+1)~
                 density,
         random=~1|cage, 
         subset(data2, thorax=="short" ))

PermTest(a1.2, B=1000)


# pos hoc for long thorax because this was significant

a1.1.1=lme(log(longevity+1)~
                 density,
         random=~1|cage, 
         subset(data2, thorax=="long" & 
                        (density=="high") |
                        (density=="isolated")))
PermTest(a1.1.1, B=1000)

a1.1.2=lme(log(longevity+1)~
                   density,
           random=~1|cage, 
           subset(data2, thorax=="long" & 
                          (density=="high") |
                          (density=="low")))
PermTest(a1.1.2, B=1000)

a1.1.3=lme(log(longevity+1)~
                   density,
           random=~1|cage, 
           subset(data2, thorax=="long" & 
                          (density=="low") |
                          (density=="isolated")))
PermTest(a1.1.3, B=1000)





