library(scales)
library(ggplot2)
getwd = getwd()
getwd
keyword = 'OTU8'
keyword = 'all'

# 
# # Assign color list to your data:
# cols <- cut(dfx$color, 2, labels = c("black", "red"))
# # ggplot() +
# #         geom_point(data = dfx,
# #                 aes(x=dfx$con_NH3, y=dfx$rm_rate,
# #                 # main= "Biomass of AOA under given NH3 concentration and remove rate",
# #                 # ylab = "remove rate", xlab = "NH3 conc",
# #                 # size = dfx$size,
# #                 # color = cols, 
# #                 # shape = 1,
# #                 pch=16))
# # + geom_point()
# # Open a pdf file
# pdf(paste("scatterplot_",keyword,".pdf",sep = ''))
# 
# Title_var = paste("Changes of ",keyword," biomass in total in the 1st iteration (hour)",sep="")
# with (dfx,(symbols(x=dfx$con_NH3, y=dfx$rm_rate,
#                    circles=dfx$size, 
#                    main= Title_var,
#                    sub = 'Under given NH3 concentration and remove rate.',
#                    xlab = "NH3 concentration (mM)", ylab = "remove rate (%)",
#                    ylim = c(0,110), 
#                    # xlim = c(0, 2.2),
#                    inches=c(0.1,1), 
#                    frame.plot = FALSE,
#                    ann=T,
#                    bg=alpha(cols, alpha = 0.3),
#                    fg=NULL
# )))
# 
# legend("topright", title = 'Change trend',
#        legend = c("reduce", "increase"),
#        pch = 19, pt.cex = 2, 
#        bty = 'o',
#        col = alpha(c('black','red'), alpha = 0.2), 
#        cex = 0.9, inset = 0.01)
# 
# legend("topleft", title = 'Weight changes (pg)',
#         legend=c("0.5", "30", "60", "90"), 
#         bty="o",
#         # inset = c(0,0.1), # where to place the legend box
#         pch=21,
#         pt.cex=c(0.5,1,2,3))
# # Close the pdf file
# dev.off() 
# 
# 
# 
# 
# pdf(paste("scatterplot_",keyword,"_ggplot.pdf",sep = ''))
# Title_var = paste("Changes of ",keyword," biomass in total in the 1st iteration (hour)",sep="")
# 
# # Change the point size, and shape
# ggplot(dfx, aes(x=con_NH3, y=rm_rate)) + 
#         geom_point(aes(color = factor(color),
#                        size = size^0.4),
#                    alpha=0.2) +
#         labs(title=Title_var, 
#              x="NH3 concentration (mM)",
#              y="remove rate (%)") +
#         scale_size(name='Weight changes'
#                    # ,labels = c("2", "0.5", "1", "0.5", "1")
#         )+
#         scale_color_manual(name='Change trend',
#                            values=c('black','red'),
#                            breaks=c("1", "2"),
#                            labels=c("Reduce", "Increase"))
# 
# dev.off() 
# 


##### Apply log transformation to x axis !! #######3
dfx = read.delim(paste(getwd,'/merge_',keyword,'_summary_for_plot.txt',sep = ''))
length(dfx$con_NH3)
length(dfx$rm_rate)
length(dfx$color)
length(dfx$size)
pdf(paste("scatterplot_",keyword,"_log_x_axis_ggplot_biomass.pdf",sep = ''))
Title_var = paste("Changes of ",keyword," biomass per hour during 10 iteration (50 hours)",sep="")

ggplot(dfx, aes(x=con_NH3, y=rm_rate)) + 
        geom_point(aes(color = factor(color),
                       size = size^0.4),
                   alpha=0.2) +
        labs(title=Title_var, 
             x="NH3 concentration (mM) with log transformation",
             y="remove rate (%)") +
        scale_size(name='Weight changes'
                   # ,labels = c("0.5", "1", "2", "3","4")
        )+
        scale_x_log10(breaks = scales::log_breaks(n = 10)) +
        annotation_logticks(sides = "b")+
        scale_color_manual(name='Change trend',
                           values=c('black','red'),
                           breaks=c("1", "2"),
                           labels=c("Reduce", "Increase")
                           # + scale_y_continuous("rm_rate", labels = as.character(rm_rate), breaks = dfx$rm_rate)
        )

dev.off() 




##### Apply log transformation to x axis !! #######3
dfx = read.delim(paste(getwd,'/merge_',keyword,'_summary_for_plot_2.txt',sep = ''))

pdf(paste("scatterplot_",keyword,"_log_x_axis_ggplot_ratio.pdf",sep = ''))
Title_var = paste("Changes of ratio of ",keyword," AOA:NOB at the end of the co-culture (the 50th hour)",sep="")

ggplot(dfx, aes(x=con_NH3, y=rm_rate)) + 
        geom_point(aes(color = factor(color),
                       size = size^0.4),
                   alpha=0.2) +
        labs(title=Title_var, 
             x="NH3 concentration (mM) with log transformation",
             y="remove rate (%)") +
        scale_size(name='Weight changes'
                   # ,labels = c("2", "0.5", "1", "0.5", "1")
        )+
        scale_x_log10(breaks = scales::log_breaks(n = 10)) +
        annotation_logticks(sides = "b")+
        scale_color_manual(name='Change trend of AOA:NOB',
                           values=c('red','black','blue'),
                           breaks=c("1", "2", "3"),
                           labels=c("Ratio no change", "Ratio reduced/increased", "NaN"))

dev.off() 


