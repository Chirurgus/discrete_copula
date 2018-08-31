# Created by Oleksandr Sorochynskyi
# On 5 6 18

library(R.matlab)

setwd("~/IdV/R")
distance_off <- readMat("data/distance.mat")$distOff
distance_on <- readMat("data/distance.mat")$distOn
distance_on_off <- readMat("data/distance.mat")$distOnOff


# cells: numbers of cells to be extracted
# binsize: size of the bins, in seconds
# cell.type: the type of cell, 1 for OFF, 6 for ON
bin_fullfield <- function(cells, bin.size = 1/60, cell.type= 1) {
  stopifnot(bin.size > 0, cells > 0);
  
  data <- readMat("data/naturalfullfield.data")
  with(data,{
  
  type <- readMat('data/twobars0703_beautiful_60Hz_2h15_30px_data.mat')$Ty[[cell.type]][[1]]
  
  stopifnot(cells <= length(type));
  
  cells <- type[cells];
  
  dt <- bin.size #inseconds
  df <-  20000*dt;
  first.rep <- 2; # Exclude the first repetion, for some reason
  n.cells = length(SpikeTimes)
  n.bins <- floor((peak.times[rep.end.index[first.rep]]-peak.times[rep.begin.index[first.rep]])/df);
  
  GoodOnes = which(rep.end.index < length(peak.times));
  rep.begin.index <- rep.begin.index[GoodOnes];
  rep.end.index <- rep.end.index[GoodOnes];
  n.reps <- length(rep.begin.index)-first.rep + 1;
  
  binrFull <- array(NA,dim= c(length(cells), n.bins, n.reps));
  for (cell in cells) {
    ST = SpikeTimes[[cell]][[1]];
    # for every repetition
    for (rep in first.rep:length(rep.begin.index)) {
      T = seq(from= peak.times[rep.begin.index[rep]],
              to= peak.times[rep.end.index[rep]],
              by= df);
      binrFull[which(cell == cells),,rep-first.rep+1] <- table(cut(ST,T));
    }
  }
  
  binrFull;
  });
}

bin_checkerboard <- function(cells, bin.size = 1/60, cell.type = 1) {
  stopifnot(bin.size > 0, cells > 0);
  
  data <- readMat('data/checkerboard.data');
  with(data, {
  
  type <- readMat('data/twobars0703_beautiful_60Hz_2h15_30px_data.mat')$Ty[[cell.type]][[1]]
  
  stopifnot(cells <= length(type));
  
  cells <- type[cells];
  
  
  df <-  20000*bin.size;
  start.offset <- 998*df;
  
  n.cells <- length(SpikeTimes);
  n.bins <- floor((peak.times[rep.end.index[1]] -
                     peak.times[rep.begin.index[1]] - start.offset)/df)
  
  
  GoodOnes <- which(rep.end.index < length(peak.times));
  
  rep.begin.index <- rep.begin.index[GoodOnes];
  rep.end.index <- rep.end.index[GoodOnes];
  n.reps <- length(rep.begin.index);
  
  binrCheck <- array(NA,dim= c(length(cells), n.bins, n.reps));
  for (cell in cells) {
    ST <- SpikeTimes[[cell]][[1]];
    # for every repetition
    for (rep in 1:n.reps) {
      T <- seq(from= peak.times[rep.begin.index[rep]]+start.offset,
              to= peak.times[rep.end.index[rep]],
              by= df);
      #T <- T[Bstart:length(T)];
      binrCheck[which(cell == cells),,rep] <- table(cut(ST,T));
    }
  }
  
  binrCheck;
  });
}


bin_barmovie <- function(cells, bar.movement= 3, bin.size= 1/60, cell.type= 1) {
  stopifnot(cells > 0,
            bar.movement > 0,
            bar.movement <= 3,
            bin.size > 0,
            cell.type > 0,
            cell.type < 8);
  data <- readMat('data/twobars0703_beautiful_60Hz_2h15_30px_data.mat') 
  
  with(data, {
  
  df <- 20000*bin.size;
  
  cells <- Ty[[cell.type]][[1]][cells]; 
  
  n.bins <- floor((Fr[rei[1,bar.movement]] - Fr[rbi[1,bar.movement]] + data$df)/df);
  n.reps <- dim(rbi)[1];
  
  binr <- array(NA, dim= c(length(cells), n.bins, n.reps));
  # For every cell
  for (cell in cells) {
    ST = Sp[[cell]][[1]];
    # for every repetition
    for (rep in 1:dim(rbi)[1]) {
      #T = Fr[rbi[rep,bar.movement]] + (0:(rei[rep,bar.movement]-rbi[rep,bar.movement]+1))*df;
      T <- seq(from= Fr[rbi[rep,bar.movement]],
               to= Fr[rei[rep,bar.movement]] + data$df,
               by= df);
      binr[which(cell == cells),,rep] <- table(cut(ST,T));
    }
  }
  binr
  });
}
