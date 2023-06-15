#capscale
#plot
#pdf(file="cap_envs_veg.pdf", width =8.5, height =6)
library(viridis)
pal<- viridis(6, option = "D") 
#plots
par(mfrow=c(2,2),mar=c(2, 2, 0.5, 0.5))
color=rgb(0,0,0,alpha=0.5) 
ordiplot(cap.env4,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env4, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env1,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env4, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)


ordiplot(cap.env1,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env1, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env2,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env1, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(cap.env2,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env2, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env3,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env2, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(cap.env3,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(cap.env3, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env4,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(cap.env3, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)



dev.off()

#plot cca
pdf(file="cca_envs_veg.pdf", width =8.5, height =6)
library(viridis)
pal<- viridis(6, option = "D") 
#plots
par(mfrow=c(2,2),mar=c(2, 2, 0.5, 0.5))
color=rgb(0,0,0,alpha=0.5) 
ordiplot(vares_cca4,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca4, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env1,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca4, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)


ordiplot(vares_cca1,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca1, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env2,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca1, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(vares_cca2,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca2, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env3,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca2, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)

ordiplot(vares_cca3,type="n")
#orditorp(cap.env,display="specie",col="red",air=0.01)
text(vares_cca3, dis="cn", scaling="sites", cex=0.5, col="red")
#orditorp(cap.env4,display="sites",cex=0.5,air=0.01)
#orditorp(cap.env,display="species",cex=0.5,air=0.01)
ordiellipse(vares_cca3, groups =metadata$Sites, col =pal, display = "sites",
            label=T, draw = "polygon", alpha = 50, label.cex=0.2, cex=0.5)



dev.off()

