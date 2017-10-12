#### "image" should be from jomovi image object
#### this function assumes that the y is called "mean", the x is "group" and the
#### tracing factor "lines" (in any)
#### errorType is an option for plotting error bars or CI. It can either  "" (empty) which
#### means no bars or CI or any other string. The string will be used as
#### the legend name and the function assumes that the data contains "lower" and "upper" varianbles
#### "theme" is passed from jmv plot function

.twoWaysPlot<-function(image,errorType,theme,depName,groupName,linesName) {
    if (errorType != '')
       dodge <- ggplot2::position_dodge(0.2)
    else
      dodge <- ggplot2::position_dodge(0)
  
    p <- ggplot2::ggplot(data=image$state$data, aes(x=group, y=mean, group=factor(lines), colour=lines)) +
       geom_line(size=.8, position=dodge) +
       labs(x=groupName, y=depName, colour=paste(linesName, errorType,sep="\n")) +
       scale_y_continuous(limits=c(min(image$state$range), max(image$state$range))) 

    if (is.factor(image$state$data$group)) {
      p <- p + geom_point(shape=21, fill='white', size=3, position=dodge)
      if (errorType != '')
          p <- p + geom_errorbar(aes(x=group, ymin=lower, ymax=upper, width=.1, group=lines), size=.8, position=dodge)
   } else {
    if (error != '')
     p <- p + geom_ribbon(aes(x=group, ymin=lower, ymax=upper,group=lines,colour=lines,fill = lines),linetype = 0,show.legend=F, alpha=.2)          
   }
 p   
}

.oneWayPlot<-function(image,errorType,theme,depName,groupName) {
  if (errorType != '')
    dodge <- ggplot2::position_dodge(0.2)
  else
    dodge <- ggplot2::position_dodge(0)
  
   p <- ggplot2::ggplot(data=image$state$data) +
        labs(x=groupName, y=depName, colour=paste("", errorType)) +
        scale_colour_manual(name=paste("", errorType), values=c(colour=theme$color[1]), labels='') +
        scale_y_continuous(limits=c(min(image$state$range), max(image$state$range))) 


  if (is.factor(image$state$data$group)) {
      p <- p + geom_point(aes(x=group, y=mean, colour='colour'), shape=21, fill=theme$fill[1], size=3)
      p <- p+geom_line(aes(x=group,y=mean,group = 1)) 
      if (errorType != '')
           p <- p + geom_errorbar(aes(x=group, ymin=lower, ymax=upper, colour='colour', width=.1), size=.8)
  }  else { 
      p <- p+geom_line(aes(x=group,y=mean)) 
      if (error != '')
      p <- p + geom_ribbon(aes(x=group, ymin=lower, ymax=upper),show.legend=F, alpha=.3)
  
  }
   p
}
