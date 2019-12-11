#################### afternoon session ####################

library(ggplot2)
library(gapminder)
library(ggrepel)
library(grid)
library(gridExtra)


file_in = '/Users/songweizhi/Desktop/gapminder_2016.csv'
plot_out = '/Users/songweizhi/Desktop/my_first_plot.png'

# read in data
gap_2016 = read.csv(file = file_in, header = T, sep = ',')


ggplot(gap_2016) +
  geom_point(aes(x = gdpPercap, y = lifeExp, color = continent)) +
  geom_smooth(aes(x = gdpPercap, y = lifeExp),
              size = 1, method = 'loess',
              color = 'red', alpha = 0.4)

gap_2016$my_labels <- gap_2016$country
gap_2016$my_labels[gap_2016$continent!= "Oceania"] <- NA

myplot = ggplot(gap_2016) +
  geom_point(aes(x=gdpPercap,y=lifeExp, color=continent, size=pop, alpha = 0.05)) +
  scale_size_continuous(range=c(0.5,20),name = "Population") + 
  labs(title = 'Gapminder DF 2016', 
       x = 'GDP/capita',
       y = 'Life Expectancy (years)') +
  theme(legend.position = 'right',
        plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(size = 10, face = 'bold')) +
  scale_color_discrete(name = 'Continent') +
  coord_cartesian(xlim = c(0, 150000), ylim =c(30, 90)) +
  scale_x_continuous(labels = c('0k', '50k', '100k', '150k')) +
  scale_y_continuous(breaks = seq(40, 85, 10)) +
  geom_text_repel(aes(x=gdpPercap, 
                      y = lifeExp, 
                      color = continent,
                      label = gap_2016$my_labels),
                  size = 3.5,
                  nudge_y = 0.3)
# save to file
ggsave(plot_out, plot = myplot, width = 12, height = 8, units = 'in')



# scale_color_grey in grey scale



