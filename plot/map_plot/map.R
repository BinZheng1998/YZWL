library(maps)
library(ggmap)
library(ggplot2)
library(readxl)
library(ggrepel)  
data <- read_excel("pig_map_metadata.xlsx")
names(data)
site<-data[,c(2:6)]
names(site)
str(site)


p <- ggplot() +

  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray70') +

  geom_point(data = site, aes(x = Longitude, y = Latitude), color = 'blue', size = 1) +

  #geom_text_repel(data = site, aes(x = Longitude, y = Latitude, label = Area), hjust = -0.5, vjust = 0.5) +

  theme_bw() +
  theme(panel.grid.minor = element_blank()) +

  scale_x_continuous(breaks=c(-180,-120, -60, 0, 60, 120, 180), expand=c(0,0),labels=c('180°','120°W','60°W','0','60°E','120°E','180°'))+
  scale_y_continuous(breaks=c(-60, -30, 0, 30, 60),expand=c(0, 0), 
                     labels=c('60°S','30°S','0','30°N','60°N'))+
  labs(x='Longitude',y='Latitude',color='Region')


print(p)

library(ggsci)


x_limits <- c(-160, 160)
y_limits <- c(-80, 80)

p <- ggplot() +

  geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group), fill = 'gray80') + 
  
  geom_point(data = site, aes(x = Longitude, y = Latitude, color = Region, size = Total_Count) ) + 

  scale_color_manual(values = c('#845EC2', '#e66e21', '#1F5897', '#FFC75F', '#C34A36','#AD5E00', '#4B4453', '#008F7A', '#a4771f', '#788a15')) +
  
  scale_size_continuous(range = c(1,3), breaks = c(10, 100, 500), labels = c("10", "100", "500")) +
  
  theme_bw() +
  theme(panel.grid = element_blank(),  
        panel.grid.minor = element_blank()) + 

  scale_color_aaas()+scale_fill_aaas()  +

  labs(x = 'Longitude', y = 'Latitude', color = 'Region', size = 'Number') +

  coord_cartesian( xlim = x_limits, ylim = y_limits)


p
