rm(list=ls())

list.of.packages <- c("sf","raster","ggplot2","ggspatial","data.table","dplyr","tidyr","parallel","bdc","taxadb","traitdataform","pbapply","tidyverse","readxl","lme4","coefplot","sjPlot","sjmisc","effects","rgdal","maptools","rgeos","terra","MuMIn","rnaturalearthdata","lsmeans","GGally","tidyterra","httr","purrr","rlist","usethis","ggpubr","leaflet","here")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

################################################################################

dir.in=here("Data") #to change accordingly to the location of the data
dir.out=here("Figs") #to change. It's the repository where the results are saved


# Load data
mydataset <- read.csv2(here(dir.in,"gen_data_final_fonseca2.csv"),
                       sep=",",dec=".",h=T) #file path in GitHub: /adaptive-potential/Data

#Data selection
## Latitude data
mydatatogo <- mydataset  %>%
    dplyr::filter(Type == "LAT",
                  shift_vel_sign == "pospos" | shift_vel_sign == "negneg", # Select only shifts in the same direction of velocity
                  SHIFT != 0, # remove non-significant shifts
                  # Nucleotide_diversity > 0 # select only GD values > 0
    ) %>% 
    mutate(
        # Climate velocity
        vel = as.numeric(velocity),
        vel_abs = abs(vel),
        vel_abs_log = log(vel_abs),
        vel_abs_log1p = log1p(vel_abs),
        
        # Shift
        SHIFT = SHIFT, # to deal with zero shift
        SHIFT_abs = abs(SHIFT),
        SHIFT_abs_log = log(SHIFT_abs),
        SHIFT_abs_log1p = log1p(SHIFT_abs),
        
        # Genetic diversity
        GD = Nucleotide_diversity, 
        GD_sqrt = sqrt(GD),
        GD_log = log(GD),
        GD_log1p = log1p(GD),
        
        # Methods
        Lat = abs(Latitude),
        Lat_band = round(Lat,0),
        ID.area = scale(ID.area),
        DUR = scale(DUR),
        LogExtent = log(Extent),
        START = scale(START),
        Param = factor(Param),
        Group = factor(Group),
        spp = factor(spp), 
        ExtentF = cut(Extent,
                      ordered_result = TRUE,
                      breaks=seq(min(Extent), max(Extent), length.out=10),
                      include.lowest=TRUE),
        NtempUnitsF = cut(NtempUnits,
                          ordered_result = TRUE,
                          breaks=seq(min(NtempUnits), max(NtempUnits), length.out=10),
                          include.lowest=TRUE)) %>%
    
    dplyr::select(
        # Genetic diversity
        GD, GD_log, GD_log1p, GD_sqrt, TajimasD,
        # Shift
        SHIFT, SHIFT_abs, SHIFT_abs_log, SHIFT_abs_log1p, 
        # SHIFT_cor, SHIFT_cor_abs, SHIFT_cor_abs_log, SHIFT_cor_raw, SHIFT_abs_log_scale,
        # Velocity
        vel, vel_abs, vel_abs_log, vel_abs_log1p, 
        trend.mean,
        # Methods + Taxonomy
        Article_ID, 
        Hemisphere,
        shift_vel_sign,
        Lat, # latitudinal position where GD was collected
        Lat_band,
        DUR, Nperiodes, LogNtempUnits, NtempUnits, Extent, LogExtent, ContinuousGrain, Quality, PrAb, ExtentF, NtempUnitsF,
        Param, Group, spp, Class, Order, Family, Genus, 
        ECO, Uncertainty_Parameter, Uncertainty_Distribution, Grain_size, Data, Article_ID,
        ID
    ) 

# transform continuous variables
cont_vars <- c(1:14,18, 20:26)

mydatatogo[,cont_vars] <- lapply(mydatatogo[,cont_vars], as.numeric)
mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

##eliminating species with uncertain obs
mydatatogo=subset(mydatatogo,spp!="Agrilus planipennis")
mydatatogo=subset(mydatatogo,spp!="Chrysodeixis eriosoma")
mydatatogo=droplevels(mydatatogo)

## Filter Classes with at least 5 species per Param
n_sps = 10

test <- mydatatogo %>%
    group_by(Class,Param) %>%
    dplyr::summarise(N_spp = length(unique(spp))) %>% # how many species per parameter?
    dplyr::filter(N_spp >= n_sps) # select classes with > n_sps per param

test <- mutate(test, Class_Param = paste(Class, Param))

mydatatogo <- mydatatogo %>%
    mutate(Class_Param = paste(Class, Param)) %>%
    filter(Class_Param %in% test$Class_Param) %>%
    select(-Class_Param)

mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

## Extra fixes
# Set the reference param level to the centroid of species obs
mydatatogo$Param <- relevel(mydatatogo$Param, ref = "O") 
################################################################################

#selection of Fonseca data
gen_div <- fread(here(dir.in,'Fonseca_etal_2023_EvoLetters_harmo.txt'))  #data from Fonseca with harmonized names

gen_div <- gen_div %>% filter(!Species %in% c("Formica-fusca_complex", "Tetramorium-caespitum_x_Tetramorium_immigrans"))

## Stats
# Gen
gen_sps <- unique(gen_div$spp_new)
match_sp=unique(mydatatogo$spp)%in%gen_sps
subset(match_sp,match_sp==T)
gen_div <- gen_div %>% dplyr::filter(spp %in% unique(mydatatogo$spp))
#biov1 <- read.csv("biov1_fixednames.csv", header = T)

## Spatial pattern
mundi <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sv")
mundi2 <- rnaturalearth::ne_countries(scale = 110, returnclass = "sv")
A116p1=crop(mundi2, extent(c(-180,180,0,90)))
A116p1$agg=1
A116p1b=aggregate(A116p1,by=A116p1$agg)
#A116p1c=terra::project(A116p1b, "+proj=moll +ellps=WGS84")
mundiX <- crop(mundi, extent(c(-179.99,179.99,-90,90)))
mundi <- crop(mundi, extent(c(-179.99,179.99,-60,90)))



# get vect gen div
vect_div <- terra::vect(gen_div,
                        geom=c("Longitude","Latitude"),crs=crs(mundi))


fgdb <- "C:/Users/brunn/ShadowDrive/CreateGeodatabaseBioShifts/Data/Study_Areas_v1/Study_Areas.gdb"
fc_list <- terra::vector_layers(fgdb)
fc_list <- fc_list[which(fc_list %in% unique(mydatatogo$ID))]
#fc_list=subset(fc_list,fc_list!="A116_P1")
fc_list="A116_P1"

for(i in 1:length(fc_list)){ cat("\r",i,"from",length(fc_list))
    tmp = terra::vect(fgdb, layer = fc_list[i])
    tmp$n_sp=length(unique(mydatatogo$spp[mydatatogo$ID==fc_list[i]]))
    tmp$n_obs=length(mydatatogo$spp[mydatatogo$ID==fc_list[i]])
    tmp$ID2=fc_list[i]
    if(i==1){
        res=tmp
    }else{
        res=rbind(res,tmp)
    }
    
}
A116p1b$ID2=res$ID2
A116p1b$REGION=res$REGION
A116p1b$Shape_Length=res$Shape_Length
A116p1b$Shape_Area=res$Shape_Area
A116p1b$n_sp=res$n_sp
A116p1b$n_obs=res$n_obs
selX=A116p1b[,c('ID2','REGION','Shape_Length','Shape_Area','n_sp','n_obs')]

fc_list <- terra::vector_layers(fgdb)
fc_list <- fc_list[which(fc_list %in% unique(mydatatogo$ID))]
fc_list=subset(fc_list,fc_list!="A116_P1")
#fc_list="A116_P1"

for(i in 1:length(fc_list)){ cat("\r",i,"from",length(fc_list))
    tmp = terra::vect(fgdb, layer = fc_list[i])
    tmp$n_sp=length(unique(mydatatogo$spp[mydatatogo$ID==fc_list[i]]))
    tmp$n_obs=length(mydatatogo$spp[mydatatogo$ID==fc_list[i]])
    tmp$ID2=fc_list[i]
    if(i==1){
        res=tmp
    }else{
        res=rbind(res,tmp)
    }
    
}
sel1=res[,c('ID2','REGION','Shape_Length','Shape_Area','n_sp','n_obs')]
sel1=rbind(sel1,selX)

# ess1=data.frame(geom(sp1))
# summary(ess1)
# ess1$x[ess1$x>=179.99]=179.99
# ess1$x[ess1$x<=(-179.99)]=(-179.99)
# unique(subset(ess1,x>=179.99)$geom)
# unique(subset(ess1,x<=(-179.99))$geom)
# geom(sp1)=ess1
nc2 <- sf::st_as_sf(sel1)
nc2 <- sf::st_transform(nc2, "+proj=moll +ellps=WGS84")



hex_area = sf::st_make_grid(nc2, cellsize=400000, square=FALSE)
grid=sf::st_as_sf(hex_area)
grid$ID=1:nrow(grid)
ig = lengths(st_intersects(hex_area, nc2)) > 0
grid$int=ig
# plot the map, and add the intersected hexagons:
plot(st_geometry(nc2))
plot(hex_area[ig], col="red", add=TRUE)
v1=sf::st_intersects(hex_area,nc2)
grid$nsp=NA
for(i in 1:nrow(grid)){
    if(grid$int[i]==T){
        grid$nsp[i]=length(unique(mydatatogo$spp[(mydatatogo$ID%in%nc2$ID2[v1[[i]]])==T]))
    }
}
grid2=subset(grid,nsp>0)
plot(grid2["nsp"])

mundi=sf::st_as_sf(mundi)
nam2=st_transform(mundi,crs=st_crs(grid))

#mundiX=sf::st_as_sf(mundiX)
#nam2=st_transform(mundiX,crs=st_crs(grid))
# library(ggplot2)
# co=st_coordinates(grid2)
# library(ggspatial)
# g1=ggplot(data = nam2) +
#   geom_sf(fill = "white") + #antiquewhite1
#   geom_sf(data = grid2, aes(fill = nsp),color="white",linewidth=0.25) +
#   scale_fill_viridis_c(alpha = .65,
#                        name = "", 
#                        guide = guide_colourbar(title.position = "left",barwidth = 0.5,barheight = 10),
#                        limits=c(1,1258),
#                        breaks=c(1,250,500,750,1000,1258)) +
#   #coord_sf(xlim = c(-180,180), ylim = c(-60,90),expand = FALSE)  +
#   #xlab("Longitude") + ylab("Latitude") +
#   ggtitle("Number of species") +#, subtitle = "(2 sites in Palm Beach County, Florida)") +
#   theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
#                                         size = 0.5), panel.background = element_rect(fill = "white"),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# g1

grid2$lognsp=log(grid2$nsp)
library(ggplot2)
co=st_coordinates(grid2)
library(ggspatial)
g2=ggplot(data = nam2) +
    geom_sf(fill = "white") + #antiquewhite1
    geom_sf(data = grid2, aes(fill = lognsp),color="white",linewidth=0.25) +
    scale_fill_viridis_c(alpha = .65,
                         name = "", 
                         guide = guide_colourbar(title.position = "left",barwidth = 0.5,barheight = 10),
                         limits=c(0,7.138),
                         breaks=log(c(1,5,10,25,50,100,200,500,1258)),
                         label=c(1,5,10,25,50,100,200,500,1258)) +
    # annotation_scale(location = "bl", width_hint = 0.25) +
    # annotation_north_arrow(location = "bl", which_north = "true", 
    #                        pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
    #                        style = north_arrow_fancy_orienteering) +
    #ggtitle("B") +
    theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                          size = 0.25), panel.background = element_rect(fill = "white"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
g2

################################################################################
#same figure for Fonseca data used in the analysis
sfgen=st_as_sf(gen_div, coords = c("Longitude","Latitude"))
st_crs(sfgen)=4326
sfgen2=st_transform(sfgen, "+proj=moll +ellps=WGS84")

v1=st_intersects(hex_area,sfgen2)
v2=lengths(v1)>0
grid$nsp_gen=NA
for(i in 1:nrow(grid)){
    if(v2[i]==T){
        grid$nsp_gen[i]=length(unique(gen_div$spp[v1[[i]]]))
    }
}
grid2=subset(grid,nsp_gen>0)
plot(grid2["nsp_gen"])
grid2$lognsp_gen=log(grid2$nsp_gen)
library(ggplot2)
co=st_coordinates(grid2)
library(ggspatial)
g3=ggplot(data = nam2) +
    geom_sf(fill = "white") + #antiquewhite1
    geom_sf(data = grid2, aes(fill = lognsp_gen),color="white",linewidth=0.25) +
    #geom_sf(data=sfgen2,color="red")+
    scale_fill_viridis_c(alpha = .65,
                         name = "Number of species (log scale)", 
                         guide = guide_colourbar(title.position = "left",barwidth = 0.5,barheight = 10),
                         limits=c(0,7.138),
                         breaks=log(c(1,5,10,25,50,100,200,500,1258)),
                         label=c(1,5,10,25,50,100,200,500,1258)) +
    # annotation_scale(location = "bl", width_hint = 0.25) +
    # annotation_north_arrow(location = "bl", which_north = "true", 
    #                        pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
    #                        style = north_arrow_fancy_orienteering) +
    #ggtitle("A") +
    theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                          size = 0.25), panel.background = element_rect(fill = "white"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.title = element_text(angle = -90), 
          legend.title.align = 0.5)+
    guides(colour = guide_colourbar(title.position = "left"),linetype="none")+
    theme(legend.key.height = unit(1.6, "cm"),
          legend.title = element_text(size = 12, angle = 90),
          legend.title.align = 0.5,
          legend.direction = "vertical")
g3

#chart formatting
p1 <- g3 + 
    theme(legend.position = "none") +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))+
    labs(tag = '(a)') +
    theme(plot.tag.position = c(0.05, 0.95))
# ggtitle("Trailing edge (TE)") +labs(tag = '(a)')+
# theme(plot.title = element_text(hjust = 0.5),
#       plot.tag.position = c(0.065, 0.975))

p2 <- g2 + 
    theme(legend.position = "none")+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))+
    labs(tag = '(b)')+
    theme(plot.tag.position = c(0.05, 0.95))
# ggtitle("Centroid (CE)") +labs(tag = '(b)')+
# theme(plot.title = element_text(hjust = 0.5))


legend <- cowplot::get_legend(g3)
gg_leg<-ggpubr::as_ggplot(legend) 
gg_leg<-gg_leg + 
    theme(plot.margin = margin(0, 0, 0, 0, "cm"))
gg_leg

design="
  13
  23
"
library(patchwork)

png(here(dir.out,"fig1.png"),unit="cm",width=16,height=13,res=300)#,width=547,height=360
(p1 + p2 + gg_leg )+
    plot_layout(design=design,widths = c(1,0.2),heights=c(0.8,0.8))
#plot_annotation(tag_levels = list(c('(a)','(b)','(c)',''))) # figure tags
dev.off()


