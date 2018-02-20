library(dplyr)
library(data.table)

pre <- read.csv('index.csv.gz')
sampled <- sample_n(pre, 1000)
sampled <- sampled[order(sampled$start_ind), ]
data.file <- read.csv('data.csv.gz')
data.df <- list()
  index.df <- list()
  current_index <- 0

for (i in 1:nrow(sampled)) {
  start <- sampled$start_ind[i]
  end <- sampled$end_ind[i]
  history <- data.file[start:end,]

  data.df <- rbindlist(list(data.df, history[,c('time', 'magnitude')]))
      index.df <- rbindlist(list(index.df, data.frame(start_ind = current_index+1, end_ind = current_index+nrow(history), id = sampled$id[i])))
current_index <- current_index + nrow(history)
}
data.df <- data.frame(data.df)
  index.df <- data.frame(index.df)
write.csv(x = data.df, file = gzfile('data1k.csv.gz'), quote = F, row.names = F)
write.csv(x = index.df, file = gzfile('index1k.csv.gz'), quote = F, row.names = F)
