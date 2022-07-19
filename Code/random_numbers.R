
data<- data.frame("Nombres y apellidos"=c("Karla Nicol Hernández Puente",
                                         "Eréndira Juanita Cano Contreras",
                                         "José Juan Blancas Vázquez"), check.names = F)
data$initials<-gsub("(?<=[A-Z])[^A-Z]+", "", data$`Nombres y apellidos`, perl = TRUE) 
data$number1<- sample(1:9, 3, replace = F)
data$number2<- sample(1:9, 3, replace = F)

data<-data %>% unite("Código", initials:number2)%>% mutate_at(c("Código"), tolower)
