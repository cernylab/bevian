
#################################

library(DT)
library(shiny)
library(shinyjs)
library(ggplot2)
library(plotly)
#library(magick)
#library(readODS)
#library(crayon)

library(foreach)
library(doParallel)

# ignore output from parallel processes
cl <- makeCluster(20)
# write output from parallel processes to stdout
#cl <- makeCluster(20, outfile="")

registerDoParallel(cl)

#Delete file browse_running if it exists
if (file.exists("browse_running")) {
  file.remove("browse_running")
}

## VALUES ##
#################################
del_data_cutoff = 100 # ignore deletes_ file if containing more than del_data_cutoff records, positions have too many errors, not worth analyzing
pixels = 224 # size of well image in pixels (integer) #TODO also split to x and y?
sizeX = 10.0 # x size of well in mm
sizeY = 10.0 # y size of well in mm
cutoff_noise <- 2.0 # only velocities larger than cutoff_noise are used
cutoff_speedup <- 10 # for finding start of movement, acceleration should be larger than cutoff_speedup
skip_speedup <- 15 # search for end of movement inside skip_speedup frames from its start
cutoff_slowdown <- 10 #
skip_slowdown <- 15 #
cutoff_bins <- 50 #

cutoff_count_movements <- 50

fps = 50.0 # frames per second
pixels_per_bin = 10 # number of pixels for grid construction (integer)
#################################

fishes <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")

#################################

## FUNCTIONS ##

#################################
#################################

dirs <- list.dirs(path=".", full.names=FALSE, recursive=FALSE)
dirs <- dirs[dirs != "RESULTS_all_csv"]

process_parallel_first <- function(dir, pixels, pixels_per_bin) {
  foreach(well=1:20, .combine=rbind, .errorhandling="remove") %dopar% {

    if (well < 10) {
      well <- paste("0", well, sep="")
    }

    files2 <- list.files(path=paste("./",dir,"/",well,"/", sep=""), pattern="*_inferred.csv$", full.names=TRUE, recursive=FALSE)
    if (length(files2) != 0) {
      dataset <- read.csv(files2, sep=",", header = TRUE)

      # transform y-axis values
      dataset$head_y <- pixels - dataset$head_y
      dataset$tail_0_y <- pixels - dataset$tail_0_y

      ## PERCENTA
      incr <- pixels_per_bin
      incr <- pixels/ceiling(pixels/incr)

      hx <- dataset$head_x   # hodnota x
      hy <- dataset$head_y   # hodnota y

      le <- length(hx) # pocet pozicii rybicky (pocet riadkov vo vstupnom file)

      for (start1 in seq(0,pixels-incr,len=ceiling(pixels/incr))) {
        for (start2 in seq(0,pixels-incr,len=ceiling(pixels/incr))) {

          end1 <- start1 + incr
          end2 <- start2 + incr

          percent <- length(which(hx > start1 & hx <= end1 & hy > start2 & hy <= end2))/le*100
          dataset$perc[hx > start1 & hx <= end1 & hy > start2 & hy <= end2] <- percent
        }
      }

      ## RYCHLOSTI A VZDIALENOSTI HLAVA-TELO
      tx <- dataset$tail_0_x   # hodnota x
      ty <- dataset$tail_0_y   # hodnota y

      diffs_list <- vector("list", le)
      lens_list <- vector("list", le)

      #cat(sprintf('\t%s: starting diffs and lens for well %s\n', date(), well), file=stdout())
      for (i in 1:le) {

        hx1 <- hx[i]
        hy1 <- hy[i]

        tx1 <- tx[i]
        ty1 <- ty[i]

        hx2 <- hx[i-1]
        hy2 <- hy[i-1]

        diff <- sqrt((hx2-hx1)**2 + (hy2-hy1)**2)
        len <- sqrt((tx1-hx1)**2 + (ty1-hy1)**2)
        diffs_list[[i]] <- diff
        lens_list[[i]] <- len

      }

      diffs_list[[1]] <- diffs_list[[2]]
      lens_list[[1]] <- lens_list[[2]]

      dataset$velocities <- unlist(diffs_list)
      dataset$lengths <- unlist(lens_list)

      ## oznacenie rybicky
      dataset$rybicka <- rep(well,nrow(dataset))
      dataset$directory <- rep(dir,nrow(dataset))

      ## accceleration + sideways shifts(ERRORs in detecing fish position)
      #cat(sprintf('\t%s: starting accels and detection of errors for well %s\n', date(), well), file=stdout())

      ## make vectors
      accel <- vector("list", le)
      accel[[1]] <- 0

      dists_tail <- vector("list", le)
      dists_head <- vector("list", le)
      dist_y_head <- vector("list", le)
      dist_y_tail <- vector("list", le)

      dists_tail[[1]] <- 0
      dists_head[[1]] <- 0
      dist_y_head[[1]] <- 0
      dist_y_tail[[1]] <- 0

      dataset$predicted_error <- rep(0,nrow(dataset))  ## NEW

      for (i in 2:le) {

        accel[[i]] <- dataset$velocities[i] - dataset$velocities[i-1]

        ## this vector
        hx1 <- hx[i-1]   ## p1   - current position of the head
        hy1 <- hy[i-1]

        tx1 <- tx[i-1]   ## p2   - current position of the tail
        ty1 <- ty[i-1]

        ## next points
        hx2 <- hx[i]  ## next p1    - next position of the head
        hy2 <- hy[i]

        tx2 <- tx[i]  ## next p2    - next position of the tail
        ty2 <- ty[i]

        ## make vectors
        # see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
        ## distance of the head from previous fish head-tail line
        dist_h1_t1_nh2_ns <- abs((tx1-hx1)*(hy1-hy2)-(hx1-hx2)*(ty1-hy1))/sqrt((tx1-hx1)**2+(ty1-hy1)**2)

        ## distance of the tail from previous fish head-tail line
        dist_h1_t1_nt2_ns <- abs((tx1-hx1)*(hy1-ty2)-(hx1-tx2)*(ty1-hy1))/sqrt((tx1-hx1)**2+(ty1-hy1)**2)

        ## squared distances between frames
        hh_dist_sq <- (hx2-hx1)**2 + (hy2-hy1)**2
        tt_dist_sq <- (tx2-tx1)**2 + (ty2-ty1)**2

        ## distance traveled along fish (i-1) axis
        h1_to_hhintersection <- sqrt(hh_dist_sq - dist_h1_t1_nh2_ns**2)
        t1_to_ttintersection <- sqrt(tt_dist_sq - dist_h1_t1_nt2_ns**2)

        ## shift to left = minus sign; shift to the right = plus sign
        ## important for the difference of distances - if fish just changed position perpendicular to the original one
        if (hx2 < hx1) {
          sign1 <- -1
        } else {
          sign1 <- 1
        }

        dist_h1_t1_nh2 <- sign1 * dist_h1_t1_nh2_ns

        ## shift to left = minus sign; shift to the right = plus sign
        ## important for the difference of distances - if fish just changed position perpendicular to the original one
        if (tx2 < tx1) {
          sign2 <- -1
        } else {
          sign2 <- 1
        }

        dist_h1_t1_nt2 <- sign2 * dist_h1_t1_nt2_ns

        dists_head[[i]] <- round(dist_h1_t1_nh2, digits=2)
        dists_tail[[i]] <- round(dist_h1_t1_nt2, digits=2)

        dist_y_head[[i]] <- round(h1_to_hhintersection, digits=2)
        dist_y_tail[[i]] <- round(t1_to_ttintersection, digits=2)

        #print(paste0("head_x: ", dist_h1_t1_nh2_ns, " tail_x: ", dist_h1_t1_nt2_ns, " head_y: ", h1_to_hhintersection, " tail_y: ", t1_to_ttintersection))
        if (h1_to_hhintersection <= 2.0 & t1_to_ttintersection <= 2.0) {    ####
          if (dist_h1_t1_nh2_ns > 3.0 & dist_h1_t1_nt2_ns > 3.0 & sign1 == sign2) {
            dataset[i, ]$predicted_error <- 1    ## ERROR
          }
        }
      }

      ## write into dataset of the fish
      dataset$heads_dists <- unlist(dists_head)
      dataset$tails_dists <- unlist(dists_tail)
      dataset$dist_y_head <- unlist(dist_y_head)
      dataset$dist_y_tail <- unlist(dist_y_tail)

      ###### skip_speedup  - to define the window
      skip_speedup <- 15    ## TOTO UZ NIEKDE BUDE
      error_end_skip <- 5   ## if no more errors occur

      #cat(sprintf('\t%s: processing detected errors for well %s\n', date(), well), file=stdout())
      ## get 1s (errors), iterate through them to get "wrong assigned areas" not just single points
      dataset_err_list <- dataset[dataset$predicted_error == 1, ]

      ## if speedups list not empty, iterate through peaks
      if (nrow(dataset_err_list) != 0) {
        ## frame numbers of the errors
        dataset_err_frames <- dataset_err_list$X       ## frames_speedups

        ## end of an error "serie"  - potential end
        err_end1 <- -1                    ## window_movement_end1 <- 0

        y <- 0

        ## iterate through error points
        #cat(sprintf('%s: well %s frames %s\n', date(), well, dataset_err_frames), file=stdout())
        for (frame in dataset_err_frames) {

          ## if frame is within a window defined by the start1 and end1 then set new boundaries to find in
          ## errors continues
          if (frame <= err_end1) {
            ## move window to find error end
            V2 <- frame + error_end_skip  # new end of window - err_end1

            err_end1 <- frame + skip_speedup
            ## start and end of the error       ### movement_real_start,movement_end,movement_NOF
            assign(paste0("error_",y), cbind(V1,V2))   ## overwritten if another peak in movement

          ## new start of movement
          } else {
            y <- y+1    ## number of all movements

            V1 <- frame  ## very first_frame

            ## so far end of the movement
            err_end1 <- frame + skip_speedup     # frame + skip    ## potential end of a movement

            V2 <- frame + error_end_skip    ###   if no more error occur set "fake" errors end
            assign(paste0("error_",y), cbind(V1, V2))

          }
        }

        #print(paste0("err1:" ,error_1))

        if (y > 0) {

          list_of_errors <- error_1
          if (y > 1) {
            for (num in 2:y) {
              list_of_errors <- rbind(list_of_errors, get(paste0("error_",num)))
            }
          }

          well_num <- as.numeric(well)
          write.table(cbind(list_of_errors),file=paste0("deletes_", dir, "_", well_num,".tsv"), sep="\t", row.names=FALSE, col.names=T, quote = F, append = FALSE)

        }
      }

      accel[[1]] <- 0
      dataset$acceleration <- unlist(accel)

   } else {
     stop(paste0("Incomplete data in ", dir, " skipping"))
   }
   # return 
    dataset
  }
}

# second part
process_parallel_second <- function(dir, datares, pixels, pixels_per_bin) {
  foreach(well2=1:20, .combine=rbind, .errorhandling="remove") %dopar% {

    if (well2 < 10) {
      well <- paste("0", well2, sep="")
    } else {
    well <- well2
    }

    dataset2 <- datares[datares$rybicka == well,]

    ## PERCENTA
    incr <- pixels_per_bin
    incr <- pixels/ceiling(pixels/incr)

    hx <- dataset2$head_x   # hodnota x
    hy <- dataset2$head_y   # hodnota y

    le <- length(hx) # pocet pozicii rybicky (pocet riadkov vo vstupnom file)

    for (start1 in seq(0,pixels-incr,len=ceiling(pixels/incr))) {
      for (start2 in seq(0,pixels-incr,len=ceiling(pixels/incr))) {

        end1 <- start1 + incr
        end2 <- start2 + incr

        percent <- round(length(which(hx > start1 & hx <= end1 & hy > start2 & hy <= end2))/le*100, digits=5)
        dataset2$perc[hx > start1 & hx <= end1 & hy > start2 & hy <= end2] <- percent
      }
    }

    ## RYCHLOSTI A VZDIALENOSTI HLAVA-TELO
    tx <- dataset2$tail_0_x   # hodnota x
    ty <- dataset2$tail_0_y   # hodnota y

    diffs_list <- vector("list", le)
    lens_list <- vector("list", le)

    for (i in 1:le) {

      hx1 <- hx[i]
      hy1 <- hy[i]

      tx1 <- tx[i]
      ty1 <- ty[i]

      hx2 <- hx[i-1]
      hy2 <- hy[i-1]

      diff <- round(sqrt((hx2-hx1)**2 + (hy2-hy1)**2), digits=2)   # distance of position of head in time = velocity
      len <- round(sqrt((tx1-hx1)**2 + (ty1-hy1)**2), digits=2)    # distance of head and tail
      diffs_list[[i]] <- diff
      lens_list[[i]] <- len

    }

    diffs_list[[1]] <- diffs_list[[2]]
    lens_list[[1]] <- lens_list[[2]]

    dataset2$velocities <- unlist(diffs_list)
    dataset2$lengths <- unlist(lens_list)

    ## oznacenie rybicky
    accel <- vector("list", le)

    for (i2 in 1:le) {

      accel[[i2]] <- dataset2$velocities[i2] - dataset2$velocities[i2-1]

    }

    accel[[1]] <- 0
    dataset2$acceleration <- unlist(accel)


    #dataresx <- dataset2
    # return
    dataset2

  }
}

# run first part parallel
for (dir in dirs) {
  #dir <- gsub(".mp4", "", filex)
  incomplete_dir = F

  if (file.exists(paste0("table_", dir, ".csv")) & file.exists(paste0("datadir_", dir, "_res.csv")) & file.exists(paste0("VELO_HIST_", dir, "_all.png")) & file.exists(paste0("VELO_", dir, "_all.png")) & file.exists(paste0("PERC_", dir, "_all.png"))) {
    next
  } else {
    print(paste0("Processing data from directory ", dir, " (runs only once, takes few minutes)"))
    print(paste0(date(), ": first part of data ..."))
    datares <- process_parallel_first(dir, pixels, pixels_per_bin)
  }

  if (typeof(datares) == "NULL") {
    print(paste0("No data found in directory ", dir, ", skipping."))
    next
  }

  ##############################################################################################################################################################
  print(paste0(date(), ": plotting figures ..."))
  ### FIGs using datares_    data
  ## PERCENTS

  png(filename = paste("PERC_", gsub("%", "%%", dir), "_all.png", sep=""),width = 3000, height = 2400)
  par(mfrow=c(4,5), oma = c(1, 1, 1, 1))

  incr <- pixels_per_bin
  incr <- pixels/ceiling(pixels/incr)

  colors <- c("lightsteelblue3", "moccasin", "lightpink", "lightpink3", "indianred1", "brown1", "red", "red3" , "darkred", "black")
  cuts <- c(0,5,7,10,15,20,35,50,60,70,100)
  bbb <- cuts[1:10]

  bins <- c("9.73913", "19.47826", "29.21739", "38.95652", "48.69565", "58.43478", "68.17391", "77.91304", "87.65217", "97.3913", "107.1304", "116.8696", "126.6087", "136.3478", "146.087", "155.8261", "165.5652", "175.3043", "185.0435", "194.7826", "204.5217", "214.2609")

  for (well in 1:20) {
    if (well < 10) {
      well <- paste("0", well, sep="")
    }

    dataset2 <- datares[datares$rybicka == well,]

    plot(c(0,pixels), c(0,pixels), xlim=c(0,pixels), ylim=c(0,pixels), col="white", axes = FALSE,ann=FALSE)
    par(new=T)
    abline(v=bins, col="lightgray", lwd=1.5)
    abline(h=bins, col="lightgray", lwd=1.5)
    par(new=T)

    for (i in 1:10) {

      hx1 <- dataset2$head_x[dataset2$perc <= cuts[i+1] & dataset2$perc > cuts[i]]
      hy1 <- dataset2$head_y[dataset2$perc <= cuts[i+1] & dataset2$perc > cuts[i]]

      plot(hx1, hy1, main = dir, xlim=c(0,pixels), ylim=c(0,pixels), col=colors[i], pch=19, cex=2.0)
      par(new=T)

    }

    abline(v=0, col="blue4", lwd=5)
    abline(v=pixels, col="blue4", lwd=5)
    abline(h=0, col="blue4", lwd=5)
    abline(h=pixels, col="blue4", lwd=5)
    par(new=F)

  }

  # umiestnenie
  par(xpd=TRUE)
  legend(220, 240, bbb, cex=1.0, col=colors, pch=21, lty=1:3)

  mtext(dir, outer=TRUE,  cex=2, line=-1.5)
  dev.off()

  ## velocities, distances head-tail
  png(filename = paste("VELO_", gsub("%", "%%", dir), "_all.png", sep=""),width = 3000, height = 2400)
  par(mfrow=c(4,5), oma = c(1, 1, 1, 1))

  for (well in 1:20) {
    if (well < 10) {
      well <- paste("0", well, sep="")
    }

    dataset3 <- datares[datares$rybicka == well,]
    velocities <- dataset3$velocities
    lengths <- dataset3$lengths

    plot(velocities, type="l", ylim=c(0,30), col="blue")
    par(new=T)
    plot(lengths, type="l", ylim=c(0,30), col="red")
  }

  par(xpd=TRUE)
  legend(220, 240, c("velo", "len"), cex=1.0, col=c("blue", "red"), pch=21, lty=1:3)

  mtext(dir, outer=TRUE,  cex=2, line=-1.5)
  dev.off()

  ## histogram - velocities
  png(filename = paste("VELO_HIST_", gsub("%", "%%", dir), "_all.png", sep=""),width = 3000, height = 2400)
  par(mfrow=c(4,5), oma = c(1, 1, 1, 1))

  for (well in 1:20) {
    if (well < 10) {
      well <- paste("0", well, sep="")
    }

    dataset4 <- datares[datares$rybicka == well,]
    velocities <- dataset4$velocities
    lengths <- dataset4$lengths

    velocities[velocities > 5.0] <- 5.0
    hist(velocities, breaks=seq(0,5.0,0.1), xlim=c(0,5))
  }

  mtext(dir, outer=TRUE,  cex=2, line=-1.5)
  dev.off()

  #### ROUNDING ####
  datares$head_x <- round(datares$head_x, digits=0)
  datares$head_y <- round(datares$head_y, digits=0)
  datares$tail_0_x <- round(datares$tail_0_x, digits=0)
  datares$tail_0_y <- round(datares$tail_0_y, digits=0)
  #### ROUNDING ####

  # run second part parallel
  print(paste0(date(), ": second part of data ..."))
  dataresx <- process_parallel_second(dir, datares, pixels, pixels_per_bin)

  ### WRITE "NEW" DATASET  ###
  write.csv(dataresx, paste("datadir_", dir, "_res.csv", sep=""), row.names = FALSE)

  ###  TABULKA  ###
  print(paste0(date(), ": populating tables ..."))

  datares1 <- dataresx[dataresx$velocities > cutoff_noise,]

  column1 <- c("ryba","means","sums","accel_counts","","ryba ","means ","sums ","accel_counts "," ","ryba  ","means  ","sums  ","accel_counts  ","  ","ryba   ","means   ","sums   ","accel_counts   ")

  for (ryb in 1:20) {
    if (ryb < 10) {
      ryb2 <- paste("0",ryb,sep="")
    } else {
      ryb2 <- ryb
    }

    datafish <- datares1[datares1$rybicka == ryb2,]   ## moze byt datares1 prazdny???  pridat podmienku
    if (nrow(datafish) != 0) {

      speedups_list <- datafish[datafish$acceleration >= cutoff_speedup,]    ## FUNKCIU

      if (nrow(speedups_list) == 0) {

        speedups_count_all <- 0
        assign(paste("accel_counts",ryb , sep=""), speedups_count_all)

      } else {

        frames_speedups <- speedups_list$X
        window_speedup_end <- 0

        for (frame in frames_speedups) {
          if (frame <= window_speedup_end) {
            next
          } else {
            window_speedup_start <- frame
            window_speedup_end <- frame+skip_speedup

            datas_in_window <- datafish[datafish$X >= window_speedup_start & datafish$X <= window_speedup_end,]
            max_speedup <- head(datas_in_window[datas_in_window$acceleration == max(datas_in_window$acceleration),]$X,1)

            if (frame == frames_speedups[1]) {          #####   ?????????????????  DOBRA PODMIENKA  ?????????????
              speedups_slowdowns <- c(max_speedup)
            } else {
              speedups_slowdowns <- cbind(speedups_slowdowns, max_speedup)
            }
          }
        }
        speedups_count_all <- length(speedups_slowdowns)
        assign(paste("accel_counts",ryb , sep=""), speedups_count_all)
      }

      means <- toString(round(mean(datafish$velocities), digits=2))
      sums <- toString(as.integer(sum(datafish$velocities)))

      assign(paste("means",ryb , sep=""), means)
      assign(paste("sums",ryb , sep=""), sums)

    } else {
      assign(paste("accel_counts",ryb , sep=""), "NA")
      assign(paste("means",ryb , sep=""), "NA")
      assign(paste("sums",ryb , sep=""), "NA")
    }
  }

  datafin <- data.frame("col1"=c("R1",means1,sums1,accel_counts1,"","R6",means6,sums6,accel_counts6,"","R11",means11,sums11,accel_counts11,"","R16",means16,sums16,accel_counts16),"col2"=c("R2",means2,sums2,accel_counts2,"","R7",means7,sums7,accel_counts7,"","R12",means12,sums12,accel_counts12,"","R17",means17,sums17,accel_counts17),"col3"=c("R3",means3,sums3,accel_counts3,"","R8",means8,sums8,accel_counts8,"","R13",means13,sums13,accel_counts13,"","R18",means18,sums18,accel_counts18),"col4"=c("R4",means4,sums4,accel_counts4,"","R9",means9,sums9,accel_counts9,"","R14",means14,sums14,accel_counts14,"","R19",means19,sums19,accel_counts19),"col5"=c("R5",means5,sums5,accel_counts5,"","R10",means10,sums10,accel_counts10,"","R15",means15,sums15,accel_counts15,"","R20",means20,sums20,accel_counts20), row.names = column1, check.rows = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
  write.csv(datafin,file=paste("table_", dir,".csv",sep=""))

  print(paste0(date(), ": directory ", dir, " DONE"))
}

# stop cluster here
stopCluster(cl)

#### LODING DATA
pathx <- "./"

####### BINS #######
incr <- pixels_per_bin
incr <- pixels/ceiling(pixels/incr)

for (start1 in seq(incr,pixels-incr,len=ceiling(pixels/incr))) {
  if (start1 == incr) {

    bins_list <- start1

  } else {

    bins_list <- cbind(bins_list, start1)

  }
}

#######################
####### CHOICES #######
print(paste0(date(), ": Loading pre-processed data (takes less than a minute) ..."))

files <- list.files(path=pathx, pattern="*_res.csv$", full.names=TRUE, recursive=FALSE)

for (filex in files) {
  name1 <- gsub("_res.csv", "", filex)
  name2 <- gsub(pathx, "", name1)
  name <- gsub("/", "", name2)
  adr <- gsub("datadir_", "", name)

  if (filex == files[1]) {
    choices_dir <- adr
  } else {
    choices_dir <- rbind(choices_dir,adr)
  }

  path <- paste(pathx, adr, "/", sep="")

  print(paste0(date(), ": reading ", filex))

  datares <- read.csv(filex, sep=",", header = TRUE)
  assign(paste("datares_", adr, sep=""), datares)

  #if (filex == files[1]) {
  #  lines_file <- length(datares[datares$rybicka == 20,]$rybicka)
  #}
  lines_file <- 10
}

print(paste0(date(), ": writing general_cutoffs.csv"))

#### ZAPISAT FILE S CUTOFFMI ####
header_col <- cbind("cutoff_name","cutoff_value")
cutoff_speedup_col <- cbind("cutoff_speedup", cutoff_speedup)
skip_speedup_col <- cbind("skip_speedup", skip_speedup)
cutoff_slowdown_col <- cbind("cutoff_slowdown", cutoff_slowdown)
skip_slowdown_col <- cbind("skip_slowdown", skip_slowdown)
cutoff_noise_col <- cbind("cutoff_noise", cutoff_noise)
cutoff_bins_col <- cbind("cutoff_bins", cutoff_bins)
condition_fps_col <- cbind("fps", fps)

table_of_cutoffs <- rbind(header_col, cutoff_speedup_col, skip_speedup_col, cutoff_slowdown_col, skip_slowdown_col, cutoff_noise_col, cutoff_bins_col, condition_fps_col)
write.csv(table_of_cutoffs, "general_cutoffs.csv")
#### ZAPISAT FILE S CUTOFFMI ####

# save a file with initial conditions (cutoffs, fish admited or not) if such file doesn't exist
# if file doesn't exist default cutoffs are set for each fish in each directory

print(paste0(date(), ": checking RESULTS_conditions.csv"))

#######  OPRAVIT !!!!!!   ####################
if (file.exists("RESULTS_conditions.csv")) {

#  "1","1h","1","10","15","10","15","1.5","ACCEPTED"
#  "2","1h","2","10","15","10","15","1.5","ACCEPTED"

#  results_cuttofs <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE)
#  head(results_cuttofs)

#  for (directory in choices_dir) {

#    dir_data <- results_cuttofs[results_cuttofs$directory == directory, ]
#    if (nrow(dir_data) == 0) {

#    print(paste0('Adding directory ', directory, ' to RESULTS_conditions.csv'))
#    cutoffs_tog <- cbind(cutoff_speedup,skip_speedup,cutoff_slowdown,skip_slowdown,cutoff_noise,cutoff_bins)

#      for (fish in fishes) {

        ## default cutoffs + if the file with selected fish (to keep) exists then use it. Fish in file keep ("ACCEPTED"), not in file "delete" ("NO")
#        if (file.exists(paste0("sel_fish_",directory,".csv"))) {
#          sel_fish_list1 <- read.csv(paste0("sel_fish_",directory,".csv"), sep=",", header=TRUE)
#          sel_fish_list <- sel_fish_list1$x

          ## all fishes in the sel_fish_file should have "ACCEPTED" in the last column
#          if (fish %in% sel_fish_list) {

#            cutoffs_fish <- as.data.frame(cbind(directory,fish,cutoffs_tog,"ACCEPTED"))

#          } else {    # fish not in the file => fish deleted => all values "NO"
#            cutoffs_fish <- as.data.frame(cbind(directory,fish,cutoffs_tog,"NO"))
#          }
#        } else {    # file with list of fishes does not exist => all 20 were accepted
#          cutoffs_fish <- as.data.frame(cbind(directory,fish,cutoffs_tog,"ACCEPTED"))
#        }
#        results_cuttofs <- as.data.frame(rbind(results_cuttofs, cutoffs_fish))
#      }

#      results_cuttofs <- as.data.frame(results_cuttofs[-1,])
#      cutoff_head <- list("directory","fish","cutoff_speedup","skip_speedup","cutoff_slowdown","skip_slowdown","cutoff_noise","cutoff_bins","acc_info")
#      names(results_cuttofs) <- cutoff_head
#      write.csv(results_cuttofs,"RESULTS_conditions.csv")

#    } else {

      print(paste0(date(), ": file with cutoffs RESULTS_conditions.csv exists"))

#    }
#  }
} else {

  print(paste0(date(), ": writing RESULTS_conditions.csv file"))

  ### cutoffs
  cutoffs_tog <- cbind(cutoff_speedup,skip_speedup,cutoff_slowdown,skip_slowdown,cutoff_noise,cutoff_bins)
  results_cuttofs <- cbind("directory","fish","cutoff_speedup","skip_speedup","cutoff_slowdown","skip_slowdown","cutoff_noise","cutoff_bins","acc_info","fps")  # cutoffs_fish

# write all cutoffs and if fish is "accepted" or not (ACCEPTED/NO)
  for (directory in choices_dir) {
    for (fish in fishes) {

      ## default cutoffs + if the file with selected fish (to keep) exists then use it. Fish in file keep ("ACCEPTED"), not in file "delete" ("NO")
      if (file.exists(paste0("sel_fish_",directory,".csv"))) {
        sel_fish_list1 <- read.csv(paste0("sel_fish_",directory,".csv"), sep=",", header=TRUE)
        sel_fish_list <- sel_fish_list1$x

        ## all fishes in the sel_fish_file should have "ACCEPTED" in the last column
        if (fish %in% sel_fish_list) {

          cutoffs_fish <- cbind(directory,fish,cutoffs_tog,"ACCEPTED",fps)

        } else {    # fish not in the file => fish deleted => all values "NO"
          cutoffs_fish <- cbind(directory,fish,cutoffs_tog,"NO",fps)
        }
      } else {    # file with list of fishes does not exist => all 20 were accepted
        cutoffs_fish <- cbind(directory,fish,cutoffs_tog,"ACCEPTED",fps)
      }
      results_cuttofs <- rbind(results_cuttofs, cutoffs_fish)
    }
  }
  results_cuttofs <- as.data.frame(results_cuttofs[-1,])
  cutoff_head <- list("directory","fish","cutoff_speedup","skip_speedup","cutoff_slowdown","skip_slowdown","cutoff_noise","cutoff_bins","acc_info","fps")
  names(results_cuttofs) <- cutoff_head
  write.csv(results_cuttofs,"RESULTS_conditions.csv")

}

#############
## FUNKCIE ##

updaterows <- function(dir, fishes, cutoff) {
  cutoff_noise <- cutoff
  table_csv <- read.csv("RESULTS_table.csv", sep=",", header=TRUE, row.names = 1)
  dirdata <- get(paste("datares_", dir, sep=""))

  cl <- makeCluster(20)
  registerDoParallel(cl)

  newrows <- foreach(fish=fishes, .combine=rbind, .export=c("get_results_table_data","MEAN_LEN_2","GET_DEL_DATA","GET_AFTER_DEL_DATA","cutoff_speedup","skip_speedup","cutoff_slowdown","skip_slowdown","cutoff_bins","cutoff_count_movements","pixels","pixels_per_bin","fps","del_data_cutoff","del_table2","get_number_of_movement_frames","get_mean_percents_selected_areas","get_maskalist_selected_areas","get_mean_per_area","get_crossings")) %dopar% {
    get_results_table_data(dir,fish,dirdata,cutoff_noise)
  }

  newrows <- as.data.frame(newrows, stringsAsFactors=F)

  for (fish in fishes) {
    rowname <- row.names(table_csv[table_csv$directory == dir & table_csv$fish == fish,])
    table_csv[toString(rowname),] <- newrows[newrows$directory == dir & newrows$fish == fish,]
  }

  stopCluster(cl)
  return(table_csv)
}

pow2 <- function(bla2) {  # get frame from on_click
  frame2 <- bla2[1,1]   # bla2 by mal byt riadok odpovedajuci vysledku z kliknutia -> prva hodnota je frame = obrazok z videa
}

pow3 <- function(bla3,dataset3) {  # get value head_x
  ds3 <- dataset3[dataset3$X == as.integer(bla3)-1,]
  head_x3 <- ds3[1,3]
}

pow4 <- function(bla4,dataset4) {  # get value head_y
  ds4 <- dataset4[dataset4$X == as.integer(bla4)-1,]
  head_y3 <- ds4[1,4]
}

pow5 <- function(bla5,dataset5) {  # get frame from dataset
  ds5 <- dataset5[dataset5$X == as.integer(bla5)-1,]
  fig_frame <- ds5[1,1]
}

get_fig_name <- function(bla6,dataset6) {  # get frame from dataset
  ds6 <- dataset6[dataset6$X == as.integer(bla6)-1,]
  fig_name2 <- ds6[1,2]
}

del_table2 <- function(dir,fish) {
  #print(paste0(date(), ": in del_table2 function"))
  if (file.exists(paste("deletes_", dir, "_", fish,".tsv",sep="")) == TRUE) {
    del_data <- data.frame("V1"=0,"V2"=0)
    tryCatch(
      { 
        del_data <- unique(read.csv(paste("deletes_", dir, "_", fish,".tsv",sep=""), sep="\t", header = T))
        if (nrow(del_data) > del_data_cutoff) {
          del_data <- data.frame("V1"=0,"V2"=0)
          print(paste("Over ", del_data_cutoff, " records in deletes_", dir, "_", fish,".tsv, data ignored",sep=""))
        } else {
          if ( ncol(del_data) != 2) {
            del_data <- data.frame("V1"=0,"V2"=0)
            print(paste("expecting two columns in deletes_", dir, "_", fish,".tsv, data ignored",sep=""))
          } else {
            if (any(names(del_data) != c("V1","V2"))) {
              del_data <- data.frame("V1"=0,"V2"=0)
              print(paste("expecting 'V1\tV2' header in deletes_", dir, "_", fish,".tsv, data ignored",sep=""))
            }
            if ( ! all(apply(del_data,2, is.numeric)) & nrow(del_data) != 0) {
              del_data <- data.frame("V1"=0,"V2"=0)
              print(paste("expecting two numeric columns in deletes_", dir, "_", fish,".tsv, data ignored",sep=""))
            }
          }
        }
      },
      warning = function(w) {print(paste0("warning in function del_table2: ",w))
      },
      error = function(e) {
        print(paste0("error in function del_table2: ",e))
        print(paste("deleting file deletes_", dir, "_", fish,".tsv, data ignored",sep=""))
        file.remove(paste0("deletes_", dir, "_", fish,".tsv"))
        del_data <- data.frame("V1"=0,"V2"=0)
      }
    )
    del_data <- data.frame(del_data[del_data$V1 > 0, ])
    return(del_data)
  }
}

GET_AFTER_DEL_DATA <- function(dir, fish, datadir) {
  cat(sprintf('%s: running GET_AFTER_DEL_DATA function', date()))
  dataset_nondel <- datadir[datadir$rybicka == fish,]
    if (file.exists(paste("deletes_", dir, "_", fish,".tsv",sep="")) == TRUE) {
      del_data <- del_table2(dir,fish)
      if (nrow(del_data) != 0 && nrow(del_data) < del_data_cutoff) {
        for (i in 1:nrow(del_data)) {

          start_del <- del_data[i,1]-1
          end_del <- del_data[i,2]-1

          dataset_nondel <- dataset_nondel[with(dataset_nondel, !(dataset_nondel$X >= start_del & dataset_nondel$X <= end_del)),]
        }
      cat(sprintf(', DONE %s\n', date()))
      return(dataset_nondel)
      } else {
        cat(sprintf(', no or too many records in deletes_ file, DONE %s\n', date()))
        return(dataset_nondel)
      }
    } else {
      cat(sprintf(', no deletes_ file present, DONE %s\n', date()))
      return(dataset_nondel)
    }
}

#### DELETED DATA
GET_DEL_DATA <- function(dir,fish, datadir) {
  cat(sprintf('%s: running GET_DEL_DATA function', date()))
  dataset_all <- datadir[datadir$rybicka == fish,]
  dataset_all_size = nrow(dataset_all)
  if (file.exists(paste("deletes_", dir, "_", fish,".tsv",sep="")) == TRUE) {
      del_data <- del_table2(dir,fish)
      if (nrow(del_data) != 0 && nrow(del_data) < del_data_cutoff) {
        for (i in 1:nrow(del_data)) {
          if ( i == 1 ) {
            if (del_data[i,1]-1 != 0) {
              start_nondel <- 0
              end_nondel <- del_data[i,1]-1
              dataset_all <- dataset_all[with(dataset_all, !(dataset_all$X >= start_nondel & dataset_all$X <= end_nondel-1)),]
            }
          } else if ( i == nrow(del_data) ) {
            if (del_data[i,2] <= dataset_all_size) {
              start_nondel <- del_data[i,2]
              end_nondel <- dataset_all_size-1
              dataset_all <- dataset_all[with(dataset_all, !(dataset_all$X >= start_nondel & dataset_all$X <= end_nondel-1)),]
            }
          } else {
            start_nondel <- del_data[i-1,2]
            end_nondel <- del_data[i,1]-1
            dataset_all <- dataset_all[with(dataset_all, !(dataset_all$X >= start_nondel & dataset_all$X <= end_nondel-1)),]
          }
        }
        dataset_del2 <- dataset_all
        cat(sprintf(', DONE %s\n', date()))
        return(dataset_del2)
      } else {
        dataset_del2 <- data.frame()
        cat(sprintf(', no or too many records in deletes_ file, DONE %s\n', date()))
        return(dataset_del2)
      }
  } else {
    dataset_all2 <- data.frame()
    cat(sprintf(', no deletes_ file present, DONE %s\n', date()))
    return(dataset_all2)
  }
}

#### DELETED DATA
## main table with summary of all_frames, sums, means; selected and not selected frames, sums, means
MEAN_LEN_2 <- function(dir,fish,datadir,cutoff_noise) {    ## datadir should be with all data, not after cutoff_noise

  #print(paste0(date(), ": in MEAN_LEN_2 function for well ", fish, ))
  cat(sprintf('%s: in MEAN_LEN_2 function for well %s from %s\n', date(), fish, dir))
  ## all dataframes
  dataset_all <- datadir[datadir$rybicka == fish, ]
  dataset_all_after_cutoff <- dataset_all[dataset_all$velocities > cutoff_noise,]

  ## before cutoff - always non-zero number of lines
  means <- toString(round(mean(dataset_all$velocities), digits=2))
  sums <- toString(as.integer(sum(dataset_all$velocities)))
  NOF <- toString(as.integer(length(dataset_all$velocities)))

  ## after cutoff - if all selected values are under cutoff then it has 0 lines
  if (nrow(dataset_all_after_cutoff) != 0) {
    means_after <- toString(round(mean(dataset_all_after_cutoff$velocities), digits=2))
    sums_after <- toString(as.integer(sum(dataset_all_after_cutoff$velocities)))
    NOF_after <- toString(as.integer(length(dataset_all_after_cutoff$velocities)))  ## number_of_frames_after
  } else {
    means_after <- 0
    sums_after <- 0
    NOF_after <- 0
  }

  table_mean_len_all <- cbind("all", means, sums, NOF, means_after, sums_after, NOF_after)

  ## not selected region/s
  dataset_nondel <- GET_AFTER_DEL_DATA(dir,fish,datadir)
  dataset_nondel_after <- dataset_nondel[dataset_nondel$velocities > cutoff_noise,]

  ## before cutoff - 0 lines would be only in case that all frames were selected
  if (nrow(dataset_nondel) != 0) {
    means_nd <- toString(round(mean(dataset_nondel$velocities), digits=2))
    sums_nd <- toString(as.integer(sum(dataset_nondel$velocities)))
    NOF_nd <- toString(as.integer(length(dataset_nondel$velocities)))
  } else {
    means_nd <- 0
    sums_nd <- 0
    NOF_nd <- 0
  }

  ## after cutoff - if all selected values are under cutoff then it has 0 lines
  if (nrow(dataset_nondel_after) != 0) {
    means_nd_after <- toString(round(mean(dataset_nondel_after$velocities), digits=2))
    sums_nd_after <- toString(as.integer(sum(dataset_nondel_after$velocities)))
    NOF_nd_after <- toString(as.integer(length(dataset_nondel_after$velocities)))
  } else {
    means_nd_after <- 0
    sums_nd_after <- 0
    NOF_nd_after <- 0
  }

  table_mean_len_nd <- cbind("all not selected", means_nd, sums_nd, NOF_nd, means_nd_after, sums_nd_after, NOF_nd_after)

  ## selected region/s
  dataset_del <- GET_DEL_DATA(dir,fish,datadir)
  dataset_del_after <- dataset_del[dataset_del$velocities > cutoff_noise,]   #####  !!!!!!!!! vacsi nez cutoff

  ## before cutoff - if no reagons selected, then it has 0 lines
  if (nrow(dataset_del) != 0) {
    means_d <- toString(round(mean(dataset_del$velocities), digits=2))
    sums_d <- toString(as.integer(sum(dataset_del$velocities)))
    NOF_d <- toString(as.integer(length(dataset_del$velocities)))
  } else {
    means_d <- 0
    sums_d <- 0
    NOF_d <- 0
  }

  ## after cutoff - if all selected values are under cutoff then it has 0 lines
  if (nrow(dataset_del_after) != 0) {
    means_d_after <- toString(round(mean(dataset_del_after$velocities), digits=2))
    sums_d_after <- toString(as.integer(sum(dataset_del_after$velocities)))
    NOF_d_after <- toString(as.integer(length(dataset_del_after$velocities)))
  } else {
    means_d_after <- 0
    sums_d_after <- 0
    NOF_d_after <- 0
  }

  table_mean_len_d <- cbind("all selected", means_d, sums_d, NOF_d, means_d_after, sums_d_after, NOF_d_after)

  table_mean_len <- rbind(table_mean_len_all, table_mean_len_nd, table_mean_len_d)

  cat(sprintf('%s: MEAN_LEN_2 DONE\n', date()))

  return(table_mean_len)
}

##########################
## table of not selected data = the rest after selection of parts of veloplot
MEAN_LEN_3 <- function(dir,fish,datadir,cutoff_noise) {    ## datadir should be with all data, not after cutoff_noise
  cat(sprintf('%s: running MEAN_LEN_3 function for well %s from %s', date(), fish, dir))
  if (file.exists(paste("deletes_", dir, "_", fish,".tsv",sep="")) == TRUE) {
    del_data <- del_table2(dir,fish)
    if (nrow(del_data) != 0) {
      dataset_all <- datadir[datadir$rybicka == fish, ]
      dataset_after_cutoff <- dataset_all[dataset_all$velocities > cutoff_noise,]

      ## dataset after cutoff is non-empty
      if (nrow(dataset_after_cutoff) != 0) {
        for (i in 1:nrow(del_data)) {

          ### add condition - of end_del is smaller than start_del
          start_del <- del_data[i,1]-1
          end_del <- del_data[i,2]-1
          start_sel <- del_data[i,1]
          end_sel <- del_data[i,2]
          dataset_del2 <- dataset_all[(dataset_all$X >= start_del & dataset_all$X <= end_del),]
          dataset_del_after <- dataset_after_cutoff[(dataset_after_cutoff$X >= start_del & dataset_after_cutoff$X <= end_del),]

          means <- toString(round(mean(dataset_del2$velocities), digits=2))
          sums <- toString(as.integer(sum(dataset_del2$velocities)))
          NOF <- toString(as.integer(length(dataset_del2$velocities)))

          means_after <- toString(round(mean(dataset_del_after$velocities), digits=2))
          sums_after <- toString(as.integer(sum(dataset_del_after$velocities)))
          NOF_after <- toString(as.integer(length(dataset_del_after$velocities)))

          table_mean_len_x <- cbind(paste0("selected_", i), means, sums, NOF, means_after, sums_after, NOF_after)

          table_tab <- cbind(table_mean_len_x, start_sel, end_sel)

          if (i == 1) {
            table_mean_len <- table_tab
          } else {
            table_mean_len <- rbind(table_mean_len, table_tab)
          }
        }
      } else {  # dataset after cutoff is empty, but dataset before cutoff is not
        for (i in 1:nrow(del_data)) {

          ### add condition - of end_del is smaller than start_del
          start_del <- del_data[i,1]-1
          end_del <- del_data[i,2]-1
          start_sel <- del_data[i,1]
          end_sel <- del_data[i,2]
          dataset_del2 <- dataset_all[(dataset_all$X >= start_del & dataset_all$X <= end_del),]

          means <- toString(round(mean(dataset_del2$velocities), digits=2))
          sums <- toString(as.integer(sum(dataset_del2$velocities)))
          NOF <- toString(as.integer(length(dataset_del2$velocities)))

          ## after cutoff
          means_after <- "NA"
          sums_after <- "NA"
          NOF_after <- "NA"

          table_mean_len_x <- cbind(paste0("selected_", i), means, sums, NOF, means_after, sums_after, NOF_after)

          table_tab <- cbind(table_mean_len_x, start_sel, end_sel)

          if (i == 1) {
            table_mean_len <- table_tab
          } else {
            table_mean_len <- rbind(table_mean_len, table_tab)
          }
        }
      }
    } else {
      table_mean_len <- cbind("","means","sums","number_of_frames","start_sel","end_sel")
    }
  } else {
    table_mean_len <- cbind("","means","sums","number_of_frames","start_sel","end_sel")
  }
  cat(sprintf(', DONE %s\n', date()))
  return(as.data.frame(table_mean_len, stringsAsFactors = FALSE))
}
##########################

####################################################
## table of selected data after selection of parts of veloplot
MEAN_LEN_4 <- function(dir,fish,datadir,cutoff_noise) {    ## datadir should be with all data, not after cutoff_noise
  cat(sprintf('%s: running MEAN_LEN_4 function for well %s from %s', date(), fish, dir))
  if (file.exists(paste("deletes_", dir, "_", fish,".tsv",sep="")) == TRUE) {
    del_data <- del_table2(dir,fish)
    if (nrow(del_data) != 0) {

      del_data <- del_data[order(del_data[,1]),]
      ## fish data with no cutoff_noise applied
      dataset_all <- datadir[datadir$rybicka == fish, ]
      ## fish data with the cutoff_noise applied
      dataset_all_after_cutoff <- dataset_all[dataset_all$velocities > cutoff_noise, ]

      le <- nrow(dataset_all)    ### uz by malo existovat - treba overit
      lend <- nrow(del_data)    ## number of lines in deletes_...csv file
      lend2 <- lend - 1

      if (del_data[1,1] != 1) {
        tog <- c(0, del_data[1,1]-1)
      } else if (del_data[1,2] != le) {
        tog <- c(del_data[1,2]-1,le-1)
      }

      if (lend > 1) {
        for (i in 1:lend2) {
          ## check if end of not-selected and start of not-selected are the same
          ## if same, 0 not-selected frames => skip
          if (del_data[i+1,1]-1 == del_data[i,2]-1) {
            next
          } else {
            first <- del_data[i,2]-1
            second <- del_data[i+1,1]-1
            tog1 <- c(first, second)
            tog <- rbind(tog, tog1)
          }
        }
      }

      tog2 <- c(del_data[lend,2]-1, le)
      tog <- rbind(tog, tog2)

      if (nrow(dataset_all_after_cutoff) != 0) {

        for (k in 1:nrow(tog)) {
          start_seq <- tog[k,1]
          end_seq <- tog[k,2]

          if (start_seq > end_seq) {
            next
          }

          if (k == 1) {
            dataset_rest <- dataset_all[(dataset_all$X >= start_seq & dataset_all$X < end_seq),]
            dataset_rest_after <- dataset_all_after_cutoff[(dataset_all_after_cutoff$X >= start_seq & dataset_all_after_cutoff$X < end_seq),]
            start_seq <- start_seq+1
          } else if (k == nrow(tog)) {
            dataset_rest <- dataset_all[(dataset_all$X > start_seq & dataset_all$X <= end_seq),]
            dataset_rest_after <- dataset_all_after_cutoff[(dataset_all_after_cutoff$X > start_seq & dataset_all_after_cutoff$X <= end_seq),]
            start_seq <- start_seq+2
          } else {
            dataset_rest <- dataset_all[(dataset_all$X > start_seq & dataset_all$X < end_seq),]
            dataset_rest_after <- dataset_all_after_cutoff[(dataset_all_after_cutoff$X > start_seq & dataset_all_after_cutoff$X < end_seq),]
            start_seq <- start_seq+2
          }

          means <- toString(round(mean(dataset_rest$velocities), digits=2))
          sums <- toString(as.integer(sum(dataset_rest$velocities)))
          NOF <- toString(as.integer(length(dataset_rest$velocities)))

          means_after <- toString(round(mean(dataset_rest_after$velocities), digits=2))
          sums_after <- toString(as.integer(sum(dataset_rest_after$velocities)))
          NOF_after <- toString(as.integer(length(dataset_rest_after$velocities)))

          table_mean_len_x <- cbind(paste0("non-sel_", k), means, sums, NOF, means_after, sums_after, NOF_after)

          table_tab <- cbind(table_mean_len_x, start_seq, end_seq)

          if (k == 1) {
            table_mean_len <- table_tab
          } else {
            table_mean_len <- rbind(table_mean_len, table_tab)
          }
        }
      } else {
        for (k in 1:nrow(tog)) {
          start_seq <- tog[k,1]
          end_seq <- tog[k,2]

          if (k == 1) {
            dataset_rest <- dataset_all[(dataset_all$X >= start_seq & dataset_all$X < end_seq),]
            dataset_rest_after <- dataset_all_after_cutoff[(dataset_all_after_cutoff$X >= start_seq & dataset_all_after_cutoff$X < end_seq),]
            start_seq <- start_seq+1
          } else if (k == nrow(tog)) {
            dataset_rest <- dataset_all[(dataset_all$X > start_seq & dataset_all$X <= end_seq),]
            dataset_rest_after <- dataset_all_after_cutoff[(dataset_all_after_cutoff$X > start_seq & dataset_all_after_cutoff$X <= end_seq),]
            start_seq <- start_seq+2
          } else {
            dataset_rest <- dataset_all[(dataset_all$X > start_seq & dataset_all$X < end_seq),]
            dataset_rest_after <- dataset_all_after_cutoff[(dataset_all_after_cutoff$X > start_seq & dataset_all_after_cutoff$X < end_seq),]
            start_seq <- start_seq+2
          }

          means <- toString(round(mean(dataset_rest$velocities), digits=2))
          sums <- toString(as.integer(sum(dataset_rest$velocities)))
          NOF <- toString(as.integer(length(dataset_rest$velocities)))

          means_after <- "NA"
          sums_after <- "NA"
          NOF_after <- "NA"

          table_mean_len_x <- cbind(paste0("non-sel_", k), means, sums, NOF, means_after, sums_after, NOF_after)

          table_tab <- cbind(table_mean_len_x, start_seq, end_seq)

          if (k == 1) {
            table_mean_len <- table_tab
          } else {
            table_mean_len <- rbind(table_mean_len, table_tab)
          }
        }
      }
    } else {
      table_mean_len <- cbind("","means","sums","number_of_frames","start_seq","end_seq")
    }
  } else {
    table_mean_len <- cbind("","means","sums","number_of_frames","start_seq","end_seq")
  }
  cat(sprintf(', DONE %s\n', date()))
  return(table_mean_len)
}
####################################################

get_data_table <- function(dir,session) {

  #print(paste0(date(), ": in get_data_table function"))
  cat(sprintf('%s: in get_data_table function for %s', date(), dir))

  table_input <- read.csv(paste("table_", dir,".csv",sep=""), sep=",", header = TRUE, stringsAsFactors = FALSE)    ## FALSE ???

  if (file.exists(paste0("sel_fish_",dir,".csv"))) {
    selected_fish1 <- read.csv(paste0("sel_fish_",dir,".csv"), sep=",", header=TRUE)
    selected_fish <- selected_fish1$x
  } else {
    selected_fish <- fishes
  }

  un_selected_fish <- fishes[!fishes %in% selected_fish]
  first <- c("1", "6", "11", "16")
  second <- c("2", "7", "12", "17")
  third <- c("3", "8", "13", "18")
  fourth <- c("4", "9", "14", "19")
  fifth <- c("5", "10", "15", "20")


  if (length(un_selected_fish) == 0) {
    cat(sprintf(', DONE %s\n', date()))
    return(table_input)
  } else {
    for (fish_to_remove in un_selected_fish) {   ## fish not in the file is the one which was unselected ~ removed
      if (fish_to_remove %in% c("1", "2", "3", "4", "5")) {
        row1 <- 2
        row2 <- 3
        row3 <- 4
      } else if (fish_to_remove %in% c("6", "7", "8", "9", "10")) {
        row1 <- 7
        row2 <- 8
        row3 <- 9
      } else if (fish_to_remove %in% c("11", "12", "13", "14", "15")) {
        row1 <- 12
        row2 <- 13
        row3 <- 14
      } else if (fish_to_remove %in% c("16", "17", "18", "19", "20")) {
        row1 <- 17
        row2 <- 18
        row3 <- 19
      }

      if (fish_to_remove %in% c("1", "6", "11", "16")) {
        first <- first[first != fish_to_remove]
        col <- 2
      } else if (fish_to_remove %in% c("2", "7", "12", "17")) {
        second <- second[second != fish_to_remove]
        col <- 3
      } else if (fish_to_remove %in% c("3", "8", "13", "18")) {
        third <- third[third != fish_to_remove]
        col <- 4
      } else if (fish_to_remove %in% c("4", "9", "14", "19")) {
        fourth <- fourth[fourth != fish_to_remove]
        col <- 5
      } else if (fish_to_remove %in% c("5", "10", "15", "20")) {
        fifth <- fifth[fifth != fish_to_remove]
        col <- 6
      }
      table_input[row1,col] <- "NA"
      table_input[row2,col] <- "NA"
      table_input[row3,col] <- "NA"
    }

    updateCheckboxGroupInput(session, "sel_fish1", label = "sel_fish", choices = c("1", "6", "11", "16"), selected = first)
    updateCheckboxGroupInput(session, "sel_fish2", label = "", choices = c("2", "7", "12", "17"), selected = second)
    updateCheckboxGroupInput(session, "sel_fish3", label = "", choices = c("3", "8", "13", "18"), selected = third)
    updateCheckboxGroupInput(session, "sel_fish4", label = "", choices = c("4", "9", "14", "19"), selected = fourth)
    updateCheckboxGroupInput(session, "sel_fish5", label = "", choices = c("5", "10", "15", "20"), selected = fifth)

    cat(sprintf(', DONE %s\n', date()))

    return(table_input)
  }
}

########## SPEEDUP_SLOWDOWNS ##########
get_speedups_and_slowdowns <- function(dir,fish,datadir,cutoff_speedup,skip_speedup,skip_slowdown,cutoff_noise) {     ###

  cat(sprintf('%s: running get_speedups_and_slowdowns function for well %s from %s', date(), fish, dir))

  # speedups_slowdowns <- cbind("max_accel_frame","min_decel_frame","max_accel_val","min_decel_val", "path_length")
  datafish <- datadir[datadir$rybicka == fish, ]     ## before applying cutoff_noise

  ## data after cutoff because all peaks (cutoff_speedup) should be bigger than selected cutoff_noise
  datadir_after_cutoff <- datafish[datafish$velocities > cutoff_noise, ]

  ## if no data after cutoff_noise return empty dataframe otherwise return list (dataframe) of all peaks (bigger or equal to cutoff_speedup)
  if (nrow(datadir_after_cutoff) != 0) {
    speedups_list <- datadir_after_cutoff[datadir_after_cutoff$acceleration >= cutoff_speedup,]
  } else {
    speedups_list <- datadir_after_cutoff # empty dataframe
  }

    ## if we got dataframe with peaks, iterate through the peaks (peaks' frame numbers)
    if (nrow(speedups_list) != 0) {

      ## get frames numbers of peaks
      frames_speedups <- speedups_list$X
      window_speedup_end <- 0

      for (frame in frames_speedups) {
        if (frame <= window_speedup_end) {

          next

        } else {

          window_speedup_start <- frame  ## first maxpeak frame
          window_speedup_end <- frame+skip_speedup ## start frame + skip to make a window to search for movement "progress"
          window_movement_start <- window_speedup_start - 1  ## "fake" start of the movement - we assume that it starts one frame before a peak

          datas_in_window <- datafish[datafish$X >= window_speedup_start & datafish$X <= window_speedup_end,]
          max_speedup <- head(datas_in_window[datas_in_window$acceleration == max(datas_in_window$acceleration),],1)   # first "global" maxpeak within an area
          max_speedup_frame <- max_speedup$X  # position of the first "global" maxpeak   !!!!! FIRST MAX MAXPEAK
          max_speedup_value <- toString(round(max_speedup$acceleration, digits=2))

          ## first lowest peak - we assume that fish started to decelerate at this point; it must be below 0; later maybe below cutoff_slowdown
          min_slowdown_start <- head(datas_in_window[datas_in_window$acceleration == min(datas_in_window$acceleration),],1)

          if (min_slowdown_start$acceleration < 0) {

            slowdown_window_start <- min_slowdown_start$X  ## first minimal slowdown = START - acceleration must be smaller than 0
            slowdown_window_end <- min_slowdown_start$X+skip_slowdown  ## END of the area

            slowdowns_window <- datafish[datafish$X <= slowdown_window_end & datafish$X >= slowdown_window_start,]
            min_slowdown_peak <- head(slowdowns_window[slowdowns_window$acceleration == min(slowdowns_window$acceleration),],1)

            min_slowdown_frame <- min_slowdown_peak$X
            min_slowdown_value <- toString(round(min_slowdown_peak$acceleration, digits=2))

            ## to calculate path length around the peak
            ## first point (frame) in our window where fish got below cutoff_noise = it stopped moving; if no such point, the last frame of the window is used instead
            first_under_cutoff <- head(slowdowns_window[slowdowns_window$velocities < cutoff_noise,],1)   ## < or <= ??
            if (nrow(first_under_cutoff) == 0) {
              window_movement_end <- slowdown_window_end
            } else {
              window_movement_end <- first_under_cutoff$X
            }

          } else {

            min_slowdown_frame <- "NA"
            min_slowdown_value <- "NA"

            ## to calculate path length around the peak   #### datas_in_window ???
            first_under_cutoff <- head(datas_in_window[datas_in_window$velocities < cutoff_noise,],1)   ## < or <= ??
            window_movement_end <- first_under_cutoff$X

          }

          around_peak_path <- datafish[datafish$X >= window_movement_start & datafish$X <= window_movement_end,]
          around_peak_path_length <- toString(round(sum(around_peak_path$velocities), digits=2))

          if (frame == frames_speedups[1]) {
            speedups_slowdowns <- cbind(max_speedup_frame, min_slowdown_frame, max_speedup_value, min_slowdown_value, around_peak_path_length, window_movement_start, window_movement_end)
          } else {
          speedups_slowdowns <- rbind(speedups_slowdowns, c(max_speedup_frame, min_slowdown_frame, max_speedup_value, min_slowdown_value, around_peak_path_length, window_movement_start, window_movement_end))
          }
        }
      }
    } else {

      max_speedup_frame <- "NA"
      min_slowdown_frame <- "NA"
      max_speedup_value <- "NA"
      min_slowdown_value <- "NA"
      around_peak_path_length <- "NA"
      window_movement_start <- "NA"
      window_movement_end <- "NA"
      speedups_slowdowns <- cbind(max_speedup_frame, min_slowdown_frame, max_speedup_value, min_slowdown_value, around_peak_path_length, window_movement_start, window_movement_end)
    }

  speedups_slowdowns <- as.data.frame(speedups_slowdowns, stringsAsFactors = FALSE)
  names(speedups_slowdowns) <- list("max_acc_fr","min_dec_fr","max_acc_val","min_dec_val", "move_len", "move_start", "move_end")
  row.names(speedups_slowdowns) <- NULL

  cat(sprintf(', DONE %s\n', date()))

  return(speedups_slowdowns)

}

########## SPEEDUP_SLOWDOWNS ##########
get_number_of_movement_frames <- function(dir,fish,datadir,cutoff_speedup,cutoff_count_movements,cutoff_noise) {     ###

  #cat(sprintf('%s: running get_number_of_movement_frames function for well %s from %s', date(), fish, dir))

  ## data for a fish
  datafish <- datadir[datadir$rybicka == fish, ]
  NOF_all <- nrow(datafish)  ## number of frames - all (around 30000 - used to be)

  ## peaks have to be above cutoff_noise; if there aren't any peaks do nothing (fill table with "NA" values - if making table)
  datadir_after_cutoff <- datafish[datafish$velocities > cutoff_noise, ]
  if (nrow(datadir_after_cutoff) != 0) {
    speedups_list <- datadir_after_cutoff[datadir_after_cutoff$acceleration >= cutoff_speedup,]
  } else {
    speedups_list <- datadir_after_cutoff # empty
  }

  ## if speedups list not empty, iterate through peaks
  if (nrow(speedups_list) != 0) {

    ## frame numbers of the peaks
    frames_speedups <- speedups_list$X

    ## start1 = frame of a last peak in a movement
    window_movement_start1 <- 0

    ## end of a movement window  - doesn't have to be a real end of whole movement - it's for searching the end frame of a movement
    window_movement_end1 <- 0

    y <- 0

    ## iterate through peaks
    for (frame in frames_speedups) {

      ## if frame is within a window defined by the start1 and end then set new boundaries to find in
      ## movement continues
      if (frame <= window_movement_end1) {

        ## move window to find movement end
        window_movement_start1 <- frame
        window_movement_end1 <- frame + cutoff_count_movements  # new end of movement

        ## movement window
        window_movement <- datafish[datafish$X > window_movement_start1 & datafish$X <= window_movement_end1,]
        movement_end <- head(window_movement[window_movement$velocities < cutoff_noise,]$X,1)

        if (length(movement_end) == 0) {
          movement_end <- window_movement_end1
        }

        ## number_of_frames
        movement_NOF <- nrow(datafish[datafish$X >= movement_real_start & datafish$X <= movement_end,])

        assign(paste0("movement_",y), c(movement_real_start,movement_end,movement_NOF))   ## overwritten if another peak in movement
      ## new start of movement
      } else {

        y <- y+1    ## number of all movements

        first_frame <- frame  ## first maxpeak frame
        window_movement_start_min <- first_frame - 10  ## "real" movement window start

        ## start of the movement ##
        window_movement_start2 <- datafish[datafish$X >= window_movement_start_min & datafish$X < first_frame,]
        movement_real_start <- tail(window_movement_start2[window_movement_start2$velocities < cutoff_noise,]$X,1)   # real start frame

        if (length(movement_real_start) == 0) {
          movement_real_start <- first_frame - 5
        }

        ## so far end of the movement
        window_movement_end1 <- frame+cutoff_count_movements # frame + skip    ## potential end of a movement
        window_movement <- datafish[datafish$X > first_frame & datafish$X < window_movement_end1,]
        movement_end <- head(window_movement[window_movement$velocities < cutoff_noise,]$X,1)   # potential end frame

        if (length(movement_end) == 0) {
          movement_end <- window_movement_end1
        }

        ## number_of_frames
        movement_NOF <- nrow(datafish[datafish$X >= movement_real_start & datafish$X <= movement_end,])

        assign(paste0("movement_",y), c(movement_real_start,movement_end,movement_NOF))
      }
    }

    list_of_movements <- movement_1

    if (y > 1) {
      for (num in 2:y) {
        list_of_movements <- rbind(list_of_movements, get(paste0("movement_",num)))
      }
    }
  } else {
    list_of_movements <- data.frame()
  }

  #cat(sprintf(', DONE %s\n', date()))

  list_of_movements2 <- as.data.frame(list_of_movements, stringsAsFactors = FALSE)
  list_of_movements2$V3
  if (length(list_of_movements2$V3) != 0) {
    movements_frames <- sum(list_of_movements2$V3)
    return(cbind(movements_frames,round(100.0*movements_frames/NOF_all,digits=2)))
  } else {
    return(cbind(0,round(0,digits=1)))
  }
}

########## MAIN TABLE - EXCLUDE FISH ##########
get_number_of_speedups <- function(datadir,fish,cutoff_speedup,skip_speedup) {

  cat(sprintf('%s: running get_number_of_speedups function for well %s', date(), fish))

  datafish <- datadir[datadir$rybicka == fish,]
  speedups_list <- datafish[datafish$acceleration >= cutoff_speedup,]    ## FUNKCIU

  if (nrow(speedups_list) == 0) {

    speedups_count_all <- 0
    assign(paste0("accel_counts",fish), speedups_count_all)

  } else {

    frames_speedups <- speedups_list$X
    window_speedup_end <- 0

    for (frame in frames_speedups) {
      if (frame <= window_speedup_end) {
        next
      } else {
        window_speedup_start <- frame
        window_speedup_end <- frame+skip_speedup

        datas_in_window <- datafish[datafish$X >= window_speedup_start & datafish$X <= window_speedup_end,]
        max_speedup <- head(datas_in_window[datas_in_window$acceleration == max(datas_in_window$acceleration),]$X,1)

        if (frame == frames_speedups[1]) {
          speedups_slowdowns <- c(max_speedup)
        } else {
          speedups_slowdowns <- cbind(speedups_slowdowns, max_speedup)
        }
      }
    }
    speedups_count_all <- length(speedups_slowdowns)

    cat(sprintf(', DONE %s\n', date()))

    return(speedups_count_all)
  }
}

### delete and return rows in delete table MEAN_LEN3
## https://github.com/stefaneng/Shiny-DeleteRowsDT/blob/master/app.R
#' Adds a row at a specified index
#'
#' @param df a data frame
#' @param row a row with the same columns as \code{df}
#' @param i the index we want to add row at.
#' @return the data frame with \code{row} added to \code{df} at index \code{i}
addRowAt <- function(df, row, i) {
  # Slow but easy to understand
  if (i > 1) {
    rbind(df[1:(i - 1), ], row, df[-(1:(i - 1)), ])
  } else {
    rbind(row, df)
  }
}

#' A column of delete buttons for each row in the data frame for the first column
#'
#' @param df data frame
#' @param id id prefix to add to each actionButton. The buttons will be id'd as id_INDEX.
#' @return A DT::datatable with escaping turned off that has the delete buttons in the first column and \code{df} in the other
deleteButtonColumn <- function(df, id, fish, dir, ...) {
  # function to create one action button as string
  f <- function(i) {
    # https://shiny.rstudio.com/articles/communicating-with-js.html
    as.character(actionButton(paste(id, i, sep="_"), label = NULL, icon = icon('trash'),
                              onclick = 'Shiny.setInputValue(\"deletePressed\",  this.id, {priority: "event"})'))
  }

  deleteCol <- unlist(lapply(seq_len(nrow(df)), f))
  write.table(cbind(df$start_sel,df$end_sel),file=paste0("deletes_", dir, "_", fish,".tsv"), sep="\t", row.names=FALSE, col.names=T, quote = F, append = FALSE)

  # Return a data table
  DT::datatable(cbind(delete = deleteCol, df),selection = 'single',
                # Need to disable escaping for html as string to work
                escape = FALSE,
                options = list(
                  # Disable sorting for the delete column
                  columnDefs = list(list(targets = 1, sortable = FALSE))
                ))

}

#' Extracts the row id number from the id string
#' @param idstr the id string formated as id_INDEX
#' @return INDEX from the id string id_INDEX
parseDeleteEvent <- function(idstr) {
  res <- as.integer(sub(".*_([0-9]+)", "\\1", idstr))
  if (! is.na(res)) res
}

###########################
##### WRITE CSV FILES #####
# PERCENTS - AREA SELECTION


get_mean_percents_selected_areas <- function(dir,fish, datadir) {

  #cat(sprintf('%s: running get_mean_percents_selected_areas function for well %s from %s', date(), fish, dir))

  ## get means for selected areas
  selected_areas_df <- read.csv("ryby_maska.csv", sep=",", header=FALSE)

  ## should be integer not float
  incrB <- pixels_per_bin
  incrS <- (pixels %% incrB)/2
  dim <- (pixels %/% incrB)+2

  perc_matrix <- matrix(ncol=dim,nrow=dim)

  col <- 0
  row <- 0

  datafish <- datadir[datadir$rybicka == fish,]
  numframes <- nrow(datafish)

  hx <- datafish$head_x   # hodnota x
  hy <- datafish$head_y   # hodnota y

  for (start1 in c(0,seq(incrS,pixels-incrS,by=incrB))) {
    col <- col+1
    for (start2 in c(0,seq(incrS,pixels-incrS,by=incrB))) {
      row <- row+1
      if (row == 1 || row == dim) {
        end2 <- start2 + incrS
      } else {
        end2 <- start2 + incrB
      }
      if (col == 1 || col == dim) {
        end1 <- start1 + incrS
      } else {
        end1 <- start1 + incrB
      }

      ## calculate percents for each number within the mask = each square
      percent_data <- datafish[which(hx > start1 & hx <= end1 & hy > start2 & hy <= end2),] ## a column in data with percent values
      perc_matrix[row,col] <- 100.0*nrow(percent_data)/numframes   ## matrix of percents

    }
    row <- 0
  }

  return(perc_matrix)

#  #cat(sprintf(', DONE %s\n', date()))

}

get_maskalist_selected_areas <- function(dir,fish, datadir) {

  #cat(sprintf('%s: running get_mean_percents_selected_areas function for well %s from %s', date(), fish, dir))

  ## get means for selected areas
  selected_areas_df <- read.csv("ryby_maska.csv", sep=",", header=FALSE)

  area_min_num <- min(selected_areas_df)
  area_max_num <- max(selected_areas_df)

  ## should be integer not float
  incrB <- pixels_per_bin
  incrS <- (pixels %% incrB)/2
  dim <- (pixels %/% incrB)+2

  datafish <- datadir[datadir$rybicka == fish,]
  numframes <- nrow(datafish)

  #datafish$maska <- rep("X",numframes)
  maskalist <- vector("list", numframes)
  # initialize to first mask value
  maskalist[sapply(maskalist, is.null)] <- area_min_num

  col <- 0
  row <- 0

  hx <- datafish$head_x   # hodnota x
  hy <- datafish$head_y   # hodnota y

  for (start1 in c(0,seq(incrS,pixels-incrS,by=incrB))) {
    col <- col+1
    for (start2 in c(0,seq(incrS,pixels-incrS,by=incrB))) {
      row <- row+1
      if (row == 1 || row == dim) {
        end2 <- start2 + incrS
      } else {
        end2 <- start2 + incrB
      }
      if (col == 1 || col == dim) {
        end1 <- start1 + incrS
      } else {
        end1 <- start1 + incrB
      }

      ## calculate percents for each number within the mask = each square
      percent_data <- datafish[which(hx > start1 & hx <= end1 & hy > start2 & hy <= end2),] ## a column in data with percent values

      ## what num is the area assigned to  +  assigned the num to each line in that area
      num_area <- selected_areas_df[row,col]

      if (length(percent_data) != 0) {
        maskalist[which(hx > start1 & hx <= end1 & hy > start2 & hy <= end2)] <- num_area
      }
    }
    row <- 0
  }

  datafish$maska <- unlist(maskalist)

  return(datafish$maska)

#  #cat(sprintf(', DONE %s\n', date()))

}

## get percents of a fish occurrence for each defined area in maska.csv
get_mean_per_area <- function(perc_matrix) {

  ## get means for selected areas
  selected_areas_df <- read.csv("ryby_maska.csv", sep=",", header=FALSE)

  area_min_num <- min(selected_areas_df)
  area_max_num <- max(selected_areas_df)

  for (area_num in area_min_num:area_max_num) {
    mean_per_area <- sum(perc_matrix[which(selected_areas_df == area_num)])
#    areas <- paste0(round(mean_per_area,digits=1), " % in area ", area_num)
    areas <- cbind(round(mean_per_area,digits=1), area_num)

    if (area_num == area_min_num) {
      areas_tog <- cbind(round(mean_per_area,digits=1), area_num)
    } else {
      areas_tog <- cbind(areas_tog, round(mean_per_area,digits=1), area_num)
    }
  }

  return(as.vector(areas_tog))

}

## get crossings from one area to another area - areas defined in maska.csv
get_crossings <- function(maska_as_row) {

  numframes <- length(maska_as_row)
  crossings2 <- c(paste0("first_second"))
  for (n in 2:numframes) {
    first <- maska_as_row[n-1]
    second <- maska_as_row[n]

    if (first != second) {
      crossings2 <- rbind(crossings2, cbind(paste0(first, "_",second)))
    }
  }

  if (!is.null(nrow(crossings2))) {
    crossings <- as.data.frame(crossings2[-1,])
    table_cross <- as.data.frame(table(crossings))

    for (i in 1:nrow(table_cross)) {
      if (i == 1) {
        list_crossings <- paste0(table_cross$crossings[1], ":", table_cross$Freq[1])
      } else {
        list_crossings <- paste0(list_crossings, ";", table_cross$crossings[i], ":", table_cross$Freq[i])
      }
    }
  } else {
    list_crossings <- "no_crossings"
  }
  return(list_crossings)
}


#for (directory in choices_dir) {

#  datadir <- read.csv(paste0("datadir_",directory,"_res.csv"), sep=",", header=TRUE)

#  for (fish in 1:20) {

##   ## make and write table with peaks
#    if (file.exists(paste0("PEAKS_",directory,"_",fish,".csv"))) {

#      next

#    } else {

#       speedups_slowdowns_data <- get_speedups_and_slowdowns(directory,fish,datadir,cutoff_speedup,skip_speedup,skip_slowdown,cutoff_noise)
#       write.csv(speedups_slowdowns_data, paste0("PEAKS_",directory,"_",fish,".csv"))

#    }


## make and write table with selected areas
#    if (file.exists(paste0("SELECTIONS_",directory,"_",fish,".csv"))) {
#      next
#    } else {
#      MEAN_LEN_3(directory,fish,datadir,cutoff_noise)
#      MEAN_LEN_4(directory,fish,datadir,cutoff_noise)

#      selections_data <- rbind(MEAN_LEN_3,MEAN_LEN_4)
#      write.csv(selections_data,paste0("SELECTIONS_",directory,"_",fish,".csv"))
#    }

#  }
#}

## make table with means, sums, number of frames - all, not-selected, selected


get_results_table_data <- function(directory,fish,datadir,cutoff_noise) {
  table_means_sums <- MEAN_LEN_2(directory,fish,datadir,cutoff_noise)

  movements_data <- get_number_of_movement_frames(directory,fish,datadir,cutoff_speedup,cutoff_count_movements,cutoff_noise)
  movements_frames <- movements_data[1]
  movements_perc <- movements_data[2]

  perc_matrix <- get_mean_percents_selected_areas(directory,fish, datadir)
  maskalist <- get_maskalist_selected_areas(directory,fish, datadir)

  areas_all <- get_mean_per_area(perc_matrix)
  area_last <- tail(areas_all,1)

  areas_names <- vector("list", area_last+1)

  for (area in 0:area_last) {
    pos <- area*2+1
    assign(paste0("area_", area), areas_all[pos])

    if (area == 0) {
      area_name <- paste0("area_", area)
      areas_percs <- get(area_name)
      areas_names[1] <- area_name

    } else {
      area_name <- paste0("area_", area)
      areas_percs <- cbind(areas_percs,get(area_name))
      areas_names[area+1] <- area_name
    }
  }

  crossings <- get_crossings(maskalist)

  results_fish <- cbind(directory,fish,table_means_sums[1,2],table_means_sums[1,3],table_means_sums[1,4],table_means_sums[1,5],table_means_sums[1,6],table_means_sums[1,7],table_means_sums[2,2],table_means_sums[2,3],table_means_sums[2,4],table_means_sums[2,5],table_means_sums[2,6],table_means_sums[2,7],table_means_sums[3,2],table_means_sums[3,3],table_means_sums[3,4],table_means_sums[3,5],table_means_sums[3,6],table_means_sums[3,7],movements_frames,movements_perc,areas_percs,crossings)
  return(results_fish)
}


if (file.exists("RESULTS_table.csv")) {

#  table_data <- read.csv("RESULTS_table.csv", sep=",", header=TRUE)

#  for (directory in choices_dir) {

#    table_dir <- table_data[table_data$directory == directory, ]

#    if (nrow(table_dir) == 0) {
#      print(paste0("Adding directory ", directory, " into RESULTS_table.csv"))

#      datadir <- read.csv(paste0("datadir_",directory,"_res.csv"), sep=",", header=TRUE)

#      for (fish in fishes) {

#        table_means_sums <- MEAN_LEN_2(directory,fish,datadir,cutoff_noise)
#        fish_means_sums <- cbind(directory,fish,table_means_sums[1,2],table_means_sums[1,3],table_means_sums[1,4],table_means_sums[1,5],table_means_sums[1,6],table_means_sums[1,7],table_means_sums[2,2],table_means_sums[2,3],table_means_sums[2,4],table_means_sums[2,5],table_means_sums[2,6],table_means_sums[2,7],table_means_sums[3,2],table_means_sums[3,3],table_means_sums[3,4],table_means_sums[3,5],table_means_sums[3,6],table_means_sums[3,7])
#        results_means_sums <- rbind(results_means_sums,fish_means_sums)

#      }
#    } else {

      print(paste0(date(), ": file RESULTS_table.csv exists"))

#    }
#  }
} else {

  print(paste0(date(), ": writing RESULTS_table.csv file"))

  for (directory in choices_dir) {
#  for (directory in c("1h")) {   ####################################

    #datadir <- read.csv(paste0("datadir_",directory,"_res.csv"), sep=",", header=TRUE)   ## staci len nacitat - uz su ulozene v premennych niekde
    datadir <- get(paste0("datares_", directory))

    for (fish in fishes) {

      results_table_data <- get_results_table_data(directory,fish,datadir,cutoff_noise)

      if (fish == 1 & directory == choices_dir[1]) {

       results_means_sums <- results_table_data

      } else {

        results_means_sums <- rbind(results_means_sums,results_table_data)

      }
    }
  }

  perc_matrix <- get_mean_percents_selected_areas(directory,fish, datadir)
  areas_all <- get_mean_per_area(perc_matrix)
  area_last <- tail(areas_all,1)

  areas_names <- vector("list", area_last+1)
  for (area in 0:area_last) {
    pos <- area*2+1
    assign(paste0("area_", area), areas_all[pos])

    if (area == 0) {
      area_name <- paste0("area_", area)
      areas_names[1] <- area_name

    } else {
      area_name <- paste0("area_", area)
      areas_names[area+1] <- area_name
    }
  }

  results_means_sums <- as.data.frame(results_means_sums,row.names=NULL)
  means_sums_head <- append(list("directory","fish","means_all","sums_all","NOF_all","means_all_after","sums_all_after","NOF_all_after","means_not_sel","sums_not_sel","NOF_not_sel","means_not_sel_after","sums_not_sel_after","NOF_not_sel_after","means_sel","sums_sel","NOF_sel","means_sel_after","sums_sel_after","NOF_sel_after","movements_frames","movements_perc"),areas_names)
  means_sums_head <- append(means_sums_head, "crossings")
  names(results_means_sums) <- means_sums_head

  write.table(results_means_sums, "RESULTS_table.csv", row.names = TRUE, sep=",")

}

bins_parallel <- function(cutoff_bins, cutoff_noise, fps, dataset, dir, fish) {
  foreach(start=seq(50,nrow(dataset)-cutoff_bins,by=cutoff_bins), .combine=rbind, .errorhandling="remove") %dopar% {

    start_in_seconds <- round(start/fps,digits=3)

    end <- start + cutoff_bins

    ## before cutoff - data in bin and their mean velocites and distances swam
    data_bin <- dataset[dataset$X+1 > start & dataset$X+1 <= end, ]

    mean_velo <- round(mean(data_bin$velocities),digits=2)
    dist_movement <- round(sum(data_bin$velocities),digits=2)     ## distance swam

    ## after cutoff
    data_bin_after_cut <- data_bin[data_bin$velocities > cutoff_noise, ]

    if (nrow(data_bin_after_cut) != 0) {

      mean_velo_after <- round(mean(data_bin_after_cut$velocities),digits=2)
      dist_movement_after <- round(sum(data_bin_after_cut$velocities),digits=2)     ## distance swam

    } else {

      mean_velo_after <- "NULL"
      dist_movement_after <- "NULL"

    }

    data_res_record <- cbind(dir, fish, start_in_seconds, start, end, mean_velo, dist_movement, mean_velo_after, dist_movement_after)
    data_res_record

  }
}

## write velocity means and sums into files ##
bins_velo_dists <- function(cutoff_bins, cutoff_noise, fps) {

  #print(paste0(date(), ": in bins_velo_dists function"))
  cat(sprintf('%s: running bins_velo_dists function', date()))

  if (dir.exists("RESULTS_all_csv") == FALSE) {
    dir.create("RESULTS_all_csv")
  }

  # ignore output from parallel processes
  cl <- makeCluster(20)
  # write output from parallel processes to stdout
  #cl <- makeCluster(20, outfile="")

  registerDoParallel(cl)

  # TODO read (fps)
  conditions <- read.csv("RESULTS_conditions.csv", header=TRUE, sep=",")
  #"","directory","fish","cutoff_speedup","skip_speedup","cutoff_slowdown","skip_slowdown","cutoff_noise","cutoff_bins","acc_info"
  #"1","12_01_21_opakovani","1","10","15","10","15","1.5","50","ACCEPTED"

  for (dir in choices_dir) {
    # TODO already read
    #datadir <- read.csv(paste0("datadir_", dir, "_res.csv"), sep=",", header=TRUE)
    datadir <- get(paste0("datares_", dir))

    for (fish in fishes) {
      # TODO
      condition <- conditions[conditions$directory == dir & conditions$fish == fish, ]
      prev_cutoff_bins <- condition$cutoff_bins
      prev_cutoff_noise <- condition$cutoff_noise
      prev_fps <- condition$fps

      if (file.exists(paste0("RESULTS_all_csv/", dir, "_", fish, "_bins.csv"))) {
        if (cutoff_bins != prev_cutoff_bins && cutoff_noise != prev_cutoff_noise && fps != prev_fps) {
          cat(sprintf(', well %s (parameters changed, recalculating)', fish))
          dataset <- datadir[datadir$rybicka == fish, ]
          dataset <- dataset[c("X", "velocities")]

          start <- 50  ## first 50 frames skipped
          end <- 0

          data_res_table <- bins_parallel(cutoff_bins, cutoff_noise, fps, dataset, dir, fish)
          write.csv(data_res_table, paste0("RESULTS_all_csv/", dir, "_", fish, "_bins.csv"))
        } else {
          cat(sprintf(', well %s (no change in parameters)', fish))
        }
      } else {
        cat(sprintf(', creating RESULTS_all_csv/%s_%s_bins.csv file',dir,fish))
        dataset <- datadir[datadir$rybicka == fish, ]
        dataset <- dataset[c("X", "velocities")]

        start <- 50  ## first 50 frames skipped
        end <- 0

        data_res_table <- bins_parallel(cutoff_bins, cutoff_noise, fps, dataset, dir, fish)
        write.csv(data_res_table, paste0("RESULTS_all_csv/", dir, "_", fish, "_bins.csv"))
      }
    }
  }

  # stop cluster here
  stopCluster(cl)

  #print(paste0(date(), ": bins_velo_dists DONE"))
  cat(sprintf(', DONE %s\n', date()))

}


bins_velo_dists(cutoff_bins, cutoff_noise, fps)

#######################################
################ SHINY ################

ui <- fluidPage(

  useShinyjs(),  # Include shinyjs

  titlePanel("Fishy fishy"),
    sidebarLayout(

      sidebarPanel(width=2,

        hidden(actionButton("WriteAllButton", "WRITE_ALL_DATA!", class = "btn-danger")),

        selectInput("datadir", "Vyber adresar:",
                   choices = choices_dir),
        selectInput("datafish", "Vyber rybicku:",
                  choices = c("all","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")),

                  numericInput("cutoff_speedup", "cutoff_speedup", value = cutoff_speedup, min = 0, max = 50, step = 0.1),
                  numericInput("skip_speedup", "skip_speedup", value = skip_speedup, min = 0, max = 50, step = 1),
                  numericInput("cutoff_slowdown", "cutoff_slowdown", value = cutoff_slowdown, min = 0, max = 50, step = 0.1),
                  numericInput("skip_slowdown", "skip_slowdown", value = skip_slowdown, min = 0, max = 50, step = 1),
                  numericInput("cutoff_noise", "cutoff_velocity_noise", value = cutoff_noise, min = 0, max = 50, step = 0.1),
                  numericInput("cutoff_count_movements", "cutoff_count_movements", value = cutoff_bins, min = 20, max = 500, step = 5),
                  actionButton("cutoffs_all_button", "Recalc dataset", class = "btn-primary"),
                  actionButton("cutoffs_alldir_button", "Recalc ALL datasets", class = "btn-danger"),
                  verbatimTextOutput("info2"),
                  verbatimTextOutput("info3"),
                  verbatimTextOutput("info4"),

                  numericInput("cutoff_bins", "Bin size (frames)", value = 50, min = 0, max = 3000, step = 50),       ## NEW
                  numericInput("fps", "frames per second", value = fps, min = 0, max = 2000, step = 50),       ## NEW
                  actionButton("cutoff_bins_button", "submit_bins", class = "btn-danger"),                      ## NEW
                  sliderInput("loop_video", "loop_video", min = -50, max = 50, value = c(-10, 10), step = 5),   ## NEW

                  plotOutput("maskaplot", width = "224px", height = "224px", ),

##
        conditionalPanel(
          condition = "input.datafish == 'all'",
          fluidRow(
            column(2, checkboxGroupInput("sel_fish1", label = "sel_fish", choices = c("1", "6", "11", "16"), selected = c("1", "6", "11", "16"))),
            column(2, checkboxGroupInput("sel_fish2", label = "", choices = c("2", "7", "12", "17"), selected = c("2", "7", "12", "17"))),
            column(2, checkboxGroupInput("sel_fish3", label = "", choices = c("3", "8", "13", "18"), selected = c("3", "8", "13", "18"))),
            column(2, checkboxGroupInput("sel_fish4", label = "", choices = c("4", "9", "14", "19"), selected = c("4", "9", "14", "19"))),
            column(2, checkboxGroupInput("sel_fish5", label = "", choices = c("5", "10", "15", "20"), selected = c("5", "10", "15", "20"))),
          ),
          fluidRow(
            column(4, actionButton("sel_fish_button", "submit_sel_fish", class = "btn-info")),
          ),
        ),
##
      ),
    mainPanel('"I get a fishy! Fishy, fishy, fishy!" -Darla', width=10,
    conditionalPanel(
        condition = "input.datafish == 'all'",
          fluidRow(
            column(6, plotOutput("img_all_perc")),
            column(6, plotOutput("img_all_velo"))
          ),
           fluidRow(
            column(6, plotOutput("img_all_hist")),
            column(6, tableOutput('TABLE'))
          ),
        ),

##################################################################################################################################
    conditionalPanel(
        condition = "input.datafish != 'all'",
      fluidRow(
        column(5, DT::dataTableOutput('SPEEDUPS_SLOWDOWNS')),
        column(7, plotlyOutput("veloPlot")),
      ),
      fluidRow(
        column(5, DT::dataTableOutput('MEAN_LEN')),
        column(4, plotlyOutput("veloPlot2", height = "200px")),
        column(3, plotOutput("distPlot", height = "250px", click = "plot_click")),
      ),
      fluidRow(
        column(4,
        div(style = "display: inline-block;vertical-align:center;",
          actionButton("left", label = "<<")),
        div(style = "display: inline-block;vertical-align:center;",
          sliderInput("frameryba", "Frame", min = 1, max = lines_file, value = 1, step = 1)),
        div(style = "display: inline-block;vertical-align:center;",
          actionButton("right", label = ">>")),
        ),
        column(1, offset = 1, actionButton("goButton", "Select!", class = "btn-success")),
        column(2, hidden(actionButton("WriteFishButton", "WRITE_FISH_DATA!", class = "btn-danger"))),
        column(4, sliderInput("delrange", "Delete", min = 1, max = lines_file, value = c(1, lines_file), step = 1)),
      ),
      fluidRow(
        uiOutput('undoUI'),
        column(6, DT::dataTableOutput('MEAN_LEN3')),
        column(6, DT::dataTableOutput('MEAN_LEN4')),
      ),
      fluidRow(
        column(width = 1, hidden(numericInput("lowx", value=0, label="lowx"))),
        column(1, hidden(numericInput("highx", value=0, label="highx"))),
        column(1, hidden(numericInput("lowy", value=0, label="lowy"))),
        column(1, hidden(numericInput("highy", value=0, label="highy"))),
      ),
    )
  )
 )
)

# nacitanie dat po vybrani datasetu (rybicky)
server <- function(input, output, session) {

  datadirInput <- reactive({ get(paste("datares_", input$datadir, sep="")) })

  output$img_all_perc <- renderImage({
  list(src = paste(pathx, "PERC_", input$datadir, "_all.png", sep=""), width = 500, height = 400)
  }, deleteFile = FALSE)

  output$img_all_velo <- renderImage({
  list(src = paste(pathx, "VELO_", input$datadir, "_all.png", sep=""), width = 500, height = 400)
  }, deleteFile = FALSE)

  output$img_all_hist <- renderImage({
  list(src = paste(pathx, "VELO_HIST_", input$datadir, "_all.png", sep=""), width = 500, height = 400)
  }, deleteFile = FALSE)

  output$TABLE <- renderTable(get_data_table(input$datadir, session))

#########################
##### DELETE BUTTON #####
  rv <- reactiveValues(
        data = NULL,
        deletedRows = NULL,
        deletedRowIndices = list()
  )

#############################
  observeEvent(input$goButton, {
##### DEL DATA UKLADANIE #####
    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv")) == TRUE) {
      del_data <- del_table2(input$datadir,input$datafish)
      if (nrow(del_data) != 0) {
        if ((input$delrange[1] <= del_data$V2) && (input$delrange[2] >= del_data$V1)) {
          print("Wrong selection! Overlapping regions")   # do info boxu!!
        } else {
          write.table(as.matrix(t(input$delrange),nrow=1, ncol=2),file=paste("deletes_", input$datadir, "_", input$datafish,".tsv",sep=""), row.names = F, col.names = F, sep="\t", quote = F, append = T)
        }
      }
    } else {
    write.table(as.matrix(t(input$delrange),nrow=1, ncol=2, dimnames = list(c(), c("from", "to"))),file=paste("deletes_", input$datadir, "_", input$datafish,".tsv",sep=""), row.names = F, col.names = T, sep="\t", quote = F)
    }

    # Clear the previous deletions
    rv$data <- as.data.frame(MEAN_LEN_3(input$datadir,input$datafish,datadirInput(),input$cutoff_noise))     ## DATADIR AFTER CUTOFF !!
    rv$deletedRows <- NULL
    rv$deletedRowIndices = list()

    output$MEAN_LEN <- DT::renderDataTable(MEAN_LEN_2(input$datadir,input$datafish,datadirInput(),input$cutoff_noise),selection = 'single')
    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,datadirInput(),input$cutoff_noise),selection = 'single')

    output$MEAN_LEN3 <- DT::renderDataTable(         ########
#      # Add the delete button column
      deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)
    )

    ## VELOPLOT ##
    output$veloPlot <- renderPlotly({
      datadir <- datadirInput()
      dataset <- datadir[datadir$rybicka == input$datafish,]
#    dataset$velocities[dataset$velocities < input$cutoff_noise] <- input$cutoff_noise  ## get data after cutoff_noise

      if (nrow(dataset) != 0) {
        dataset_after_del <- GET_AFTER_DEL_DATA(input$datadir,input$datafish,datadir)
        dataset_after_del$velocities[dataset_after_del$velocities < input$cutoff_noise] <- input$cutoff_noise  ## get data after cutoff_noise
        dataset_deleted <- GET_DEL_DATA(input$datadir,input$datafish,datadir)

        if (nrow(dataset_deleted) == 0) {
          dataset_deleted <- dataset[(dataset$X >= 1 & dataset$X <= 3),]
        } else {
          dataset_deleted$velocities[dataset_deleted$velocities < input$cutoff_noise] <- input$cutoff_noise  ## get data after cutoff_noise
        }

        fig <- plot_ly(dataset_deleted, x=~X, y=~lengths, type="scatter", mode="lines", line = list(color = "gray"), name="h-t dist del") %>%
        add_trace(y=~velocities, mode="lines", line = list(color = "gray"), name="velocity del") %>%
        add_trace(y=~acceleration, mode="lines", line = list(color = "gray"), name="acceleration del") %>%
        add_trace(data=dataset_after_del, y=~lengths, mode="lines", line = list(color = "red"), name="h-t dist") %>%
        add_trace(data=dataset_after_del, y=~velocities, mode="lines", line = list(color = "blue"), name="velocity") %>%
        add_trace(data=dataset_after_del, y=~acceleration, mode="lines", line = list(color = "magenta"), name="acceleration") %>%
        layout(showlegend = TRUE,  yaxis = list(hoverformat = '.1f'), xaxis = list(hoverformat = 'd')) %>%
        layout(xaxis = list(range=c(1,lines_file))) %>%
        event_register('plotly_relayout') %>%
        event_register('plotly_click')

      } else {
        faked <- rep(1,lines_file)
        fig <- plot_ly(dataset, x=1:lines_file, y=faked, type="scatter", mode="lines", line = list(color = "gray"), name="h-t dist del") %>%
        layout(xaxis = list(range=c(1,lines_file))) %>%
        layout(yaxis = list(range=c(-10,50)))
      }

    })
  ## VELOPLOT ##
  ##############

    ## update the RESULTS_table.csv file based on the cutoff_noise change
    ## the MEAN_LEN_2 function takes into account selected and not selected parts from deletes_...csv files
    ## UROBIT FUNKCIU !!!!!!!!!!!!!
    table_csv <- read.csv("RESULTS_table.csv", sep=",", header=TRUE, row.names = 1)
    rowname <- row.names(table_csv[table_csv$directory == input$datadir & table_csv$fish == input$datafish,])
    table_csv[toString(rowname),] <- get_results_table_data(input$datadir,input$datafish,datadirInput(),input$cutoff_noise)
    write.csv <- write.csv(table_csv,"RESULTS_table.csv")

  })

  observeEvent(input$deletePressed, {
    rowNum <- parseDeleteEvent(input$deletePressed)
    dataRow <- rv$data[rowNum,]

    # Put the deleted row into a data frame so we can undo
    # Last item deleted is in position 1
    rv$deletedRows <- rbind(dataRow, rv$deletedRows)
    rv$deletedRowIndices <- append(rv$deletedRowIndices, rowNum, after = 0)

    # Delete the row from the data frame
    rv$data <- rv$data[-rowNum,]
  })

  observeEvent(input$undo, {
    if(nrow(rv$deletedRows) > 0) {
      row <- rv$deletedRows[1, ]
      rv$data <- addRowAt(rv$data, row, rv$deletedRowIndices[[1]])

      # Remove row
      rv$deletedRows <- rv$deletedRows[-1,]
      # Remove index
      rv$deletedRowIndices <- rv$deletedRowIndices[-1]
    }
  })

  # Disable the undo button if we have not deleted anything
  output$undoUI <- renderUI({
    if(!is.null(rv$deletedRows) && nrow(rv$deletedRows) > 0) {
      actionButton('undo', label = 'Undo delete', icon('undo'))
    } else {
      actionButton('undo', label = 'Undo delete', icon('undo'), disabled = TRUE)
    }
  })

##### DELETE BUTTON #####
#########################

################################################3

  observeEvent(input$sel_fish_button, {
    selected_fish <- unique(c(input$sel_fish1, input$sel_fish2, input$sel_fish3, input$sel_fish4, input$sel_fish5))
    write.csv(selected_fish, paste0("sel_fish_",input$datadir,".csv"))

    if (length(selected_fish) == 0) {
      print("No fish selected! Select fish")
    } else {
      output$TABLE <- renderTable(get_data_table(input$datadir, session))
    }

    ## sel_fish
    ## change ACCEPTED to NOT_ACCEPTED
    un_selected_fish <- fishes[!fishes %in% selected_fish]
    cutoffs_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=T, row.names = 1)
    for (un_sel_fish in un_selected_fish) {
      cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == input$datafish,]$acc_info <- "NO"
      write.csv(cutoffs_csv,"RESULTS_conditions.csv")
    }
  })

###### bins_cutoff
  observeEvent(input$cutoff_bins_button, {
    bins_velo_dists(input$cutoff_bins, input$cutoff_noise, input$fps)
  })

##########################################
##### NOISE and other cutoffs BUTTON #####
  observeEvent(input$cutoffs_button, {

    #### ZAPISAT FILE S CUTOFFMI ####
    header <- cbind("cutoff_name","cutoff_value")
    cutoff_speedup_col <- cbind("cutoff_speedup", input$cutoff_speedup)
    skip_speedup_col <- cbind("skip_speedup", input$skip_speedup)
    cutoff_slowdown_col <- cbind("cutoff_slowdown", input$cutoff_slowdown)
    skip_slowdown_col <- cbind("skip_slowdown", input$skip_slowdown)
    cutoff_noise_col <- cbind("cutoff_noise", input$cutoff_noise)

    table_of_cutoffs <- rbind(header, cutoff_speedup_col, skip_speedup_col, cutoff_slowdown_col, skip_slowdown_col, cutoff_noise_col)
    write.csv(table_of_cutoffs, "general_cutoffs.csv")
    #### ZAPISAT FILE S CUTOFFMI ####

    datadir_orig <- datadirInput()
    datadir_after_cutoff <- datadir_orig[datadir_orig$velocities > input$cutoff_noise, ]

    output$MEAN_LEN <- DT::renderDataTable(MEAN_LEN_2(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')

    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv"))) {
      data <- del_table2(input$datadir,input$datafish)
      if (nrow(data) != 0) {
        output$MEAN_LEN3 <- DT::renderDataTable(
#       # Add the delete button column
        deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)
        )
      } else {
        output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
      }
    } else {
      output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
    }

    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
    output$SPEEDUPS_SLOWDOWNS <- DT::renderDataTable(get_speedups_and_slowdowns(input$datadir,input$datafish,datadir_orig,input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise),selection = 'single')

  })

##############################
##### cutoffs_all_button #####
observeEvent(input$cutoffs_all_button, {

  for (fish in fishes) {
    cutoffs_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)
    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == fish,]$cutoff_speedup <- input$cutoff_speedup
    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == fish,]$skip_speedup <- input$skip_speedup
    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == fish,]$cutoff_slowdown <- input$cutoff_slowdown
    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == fish,]$skip_slowdown <- input$skip_slowdown
    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == fish,]$cutoff_noise <- input$cutoff_noise
    write.csv <- write.csv(cutoffs_csv,"RESULTS_conditions.csv")
  }
})

  observeEvent(input$cutoffs_alldir_button, {
    for (dir in choices_dir) {
      for (fish in fishes) {
        cutoffs_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)
        cutoffs_csv[cutoffs_csv$directory == dir & cutoffs_csv$fish == fish,]$cutoff_speedup <- input$cutoff_speedup
        cutoffs_csv[cutoffs_csv$directory == dir & cutoffs_csv$fish == fish,]$skip_speedup <- input$skip_speedup
        cutoffs_csv[cutoffs_csv$directory == dir & cutoffs_csv$fish == fish,]$cutoff_slowdown <- input$cutoff_slowdown
        cutoffs_csv[cutoffs_csv$directory == dir & cutoffs_csv$fish == fish,]$skip_slowdown <- input$skip_slowdown
        cutoffs_csv[cutoffs_csv$directory == dir & cutoffs_csv$fish == fish,]$cutoff_noise <- input$cutoff_noise
        write.csv(cutoffs_csv,"RESULTS_conditions.csv")
      }
    }
  })

#########################################
# observeEvents  - cutoffs
#"","directory","fish","cutoff_speedup","skip_speedup","cutoff_slowdown","skip_slowdown","cutoff_noise","acc_info"
#"1","1h","1","10","15","10","15","1.5","ACCEPTED"
#"2","1h","2","10","15","10","15","1.5","ACCEPTED"

## cutoff_speedup
observeEvent(input$cutoff_speedup, {
  if (input$datafish != "all") {
    cutoffs_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)
    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == input$datafish,]$cutoff_speedup <- input$cutoff_speedup
    write.csv <- write.csv(cutoffs_csv,"RESULTS_conditions.csv")

    datadir_orig <- datadirInput()

    output$SPEEDUPS_SLOWDOWNS <- DT::renderDataTable(get_speedups_and_slowdowns(input$datadir,input$datafish,datadirInput(),input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise),selection = 'single')

    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv"))) {
      data <- del_table2(input$datadir,input$datafish)
      if (nrow(data) != 0) {
        output$MEAN_LEN3 <- DT::renderDataTable(
#       # Add the delete button column
        deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)
        )
      } else {
        output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
      }
    } else {
      output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
    }

    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,datadirInput(),input$cutoff_noise),selection = 'single')
  }
})


## skip_speedup
observeEvent(input$skip_speedup, {
  if (input$datafish != "all") {
    cutoffs_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)

    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == input$datafish,]$skip_speedup <- input$skip_speedup
    write.csv <- write.csv(cutoffs_csv,"RESULTS_conditions.csv")

    datadir_orig <- datadirInput()

    output$SPEEDUPS_SLOWDOWNS <- DT::renderDataTable(get_speedups_and_slowdowns(input$datadir,input$datafish,datadirInput(),input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise),selection = 'single')

    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv"))) {
      data <- del_table2(input$datadir,input$datafish)
      if (nrow(data) != 0) {
        output$MEAN_LEN3 <- DT::renderDataTable(
#       # Add the delete button column
        deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)
        )
      } else {
        output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
      }
    } else {
      output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
    }

    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,datadirInput(),input$cutoff_noise),selection = 'single')
  }
})


## cutoff_slowdown
observeEvent(input$cutoff_slowdown, {
  if (input$datafish != "all") {
    cutoffs_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)

    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == input$datafish,]$cutoff_slowdown <- input$cutoff_slowdown
    write.csv <- write.csv(cutoffs_csv,"RESULTS_conditions.csv")

    datadir_orig <- datadirInput()

    output$SPEEDUPS_SLOWDOWNS <- DT::renderDataTable(get_speedups_and_slowdowns(input$datadir,input$datafish,datadirInput(),input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise),selection = 'single')

    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv"))) {
      data <- del_table2(input$datadir,input$datafish)
      if (nrow(data) != 0) {
        output$MEAN_LEN3 <- DT::renderDataTable(
#       # Add the delete button column
        deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)    #### !!!!!!!!!!!!!!!
        )
      } else {
        output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
      }
    } else {
      output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
    }

    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,datadirInput(),input$cutoff_noise),selection = 'single')
  }
})

## skip_slowdown
observeEvent(input$skip_slowdown, {
  if (input$datafish != "all") {
    cutoffs_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)

    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == input$datafish,]$skip_slowdown <- input$skip_slowdown
    write.csv <- write.csv(cutoffs_csv,"RESULTS_conditions.csv")

    datadir_orig <- datadirInput()

    output$SPEEDUPS_SLOWDOWNS <- DT::renderDataTable(get_speedups_and_slowdowns(input$datadir,input$datafish,datadirInput(),input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise),selection = 'single')

    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv"))) {
      data <- del_table2(input$datadir,input$datafish)
      if (nrow(data) != 0) {
        output$MEAN_LEN3 <- DT::renderDataTable(
#       # Add the delete button column
        deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)
        )
      } else {
        output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
      }
    } else {
      output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
    }

    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,datadirInput(),input$cutoff_noise),selection = 'single')
  }
})

observeEvent(input$cutoffs_all_button, {
  # for all fish in only one directory - the active one
  cat(sprintf('%s: Recalculating data for ALL fish in current dataset ... ', date()))
  showModal(modalDialog("Recalculating data for ALL fish in current dataset ...", footer=NULL))
  newrows <- updaterows(input$datadir, fishes, input$cutoff_noise)
  write.csv <- write.csv(newrows,"RESULTS_table.csv")
  removeModal()
  cat(sprintf('%s: all DONE\n', date()))
})

observeEvent(input$cutoffs_alldir_button, {
  # for all fish in all directories
  cat(sprintf('%s: Recalculating data for ALL fish in ALL datasets ... ', date()))
  showModal(modalDialog("Recalculating data for ALL fish in ALL datasets ...", footer=NULL))
  for (directory in dirs) {
    newrows <- updaterows(directory, fishes, input$cutoff_noise)
    write.csv <- write.csv(newrows,"RESULTS_table.csv")
  }
  removeModal()
  cat(sprintf('%s: all DONE\n', date()))
})

## cutoff_noise
observeEvent(input$cutoff_noise, {
  if (input$datafish != "all") {
    cutoffs_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)

    cutoffs_csv[cutoffs_csv$directory == input$datadir & cutoffs_csv$fish == input$datafish,]$cutoff_noise <- input$cutoff_noise
    write.csv <- write.csv(cutoffs_csv,"RESULTS_conditions.csv")

    ## update the RESULTS_table.csv file based on the cutoff_noise change
    ## the MEAN_LEN_2 function takes into account selected and not selected parts from deletes_...csv files
    table_csv <- read.csv("RESULTS_table.csv", sep=",", header=TRUE, row.names = 1)
    rowname <- row.names(table_csv[table_csv$directory == input$datadir & table_csv$fish == input$datafish,])
    table_csv[toString(rowname),] <- get_results_table_data(input$datadir,input$datafish,datadirInput(),input$cutoff_noise)
    write.csv <- write.csv(table_csv,"RESULTS_table.csv")

    datadir_orig <- datadirInput()

    rv$data <- as.data.frame(MEAN_LEN_3(input$datadir,input$datafish,datadirInput(),input$cutoff_noise))     ## DATADIR AFTER CUTOFF !!
    rv$deletedRows <- NULL
    rv$deletedRowIndices = list()

    ## update table of speedups and slowdowns based on the cutoff_noise change
    output$SPEEDUPS_SLOWDOWNS <- DT::renderDataTable(get_speedups_and_slowdowns(input$datadir,input$datafish,datadirInput(),input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise),selection = 'single')

    ## update table with selected areas (number of frames, means, sums) based on the cutoff_noise change
    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv"))) {
      data <- del_table2(input$datadir,input$datafish)
      if (nrow(data) != 0) {

        output$MEAN_LEN3 <- DT::renderDataTable(
#       # Add the delete button column if files with selected areas exists and is not empty  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)
        )
      } else {
        output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
      }
    } else {
      output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')

    }

    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,datadirInput(),input$cutoff_noise),selection = 'single')

  }
})

##### NOISE BUTTON #####
########################

## speedup selected from table
observeEvent(input$SPEEDUPS_SLOWDOWNS_rows_selected, {

  datadir_orig <- datadirInput()

  #tablex <- as.data.frame(get_speedups_and_slowdowns(input$datadir,input$datafish,datadir_orig,input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise), stringsAsFactors = FALSE)
  tablex <- get_speedups_and_slowdowns(input$datadir,input$datafish,datadir_orig,input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise)
  #print(tablex[input$SPEEDUPS_SLOWDOWNS_rows_selected,])
  peak_position <- as.integer(tablex[input$SPEEDUPS_SLOWDOWNS_rows_selected,1])

  updateSliderInput(session, "frameryba", value = peak_position)

})

## error selected from deletes_ table
observeEvent(input$MEAN_LEN3_rows_selected, {

  datadir_orig <- datadirInput()

  tablex <- MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise)

  selected_row <- tablex[input$MEAN_LEN3_rows_selected, ]

  start_position <- as.integer(selected_row$start_sel)
  end_position <- as.integer(selected_row$end_sel)

  updateSliderInput(session, "frameryba", value = start_position)

  updateSliderInput(session, "delrange", value = c(start_position,end_position),  min = start_position, max = end_position, step = 1)

})

###########################################################################################################################################
observeEvent(input$left, {
    updateSliderInput(session, "frameryba", value = input$frameryba - 10)
})
observeEvent(input$right, {
    updateSliderInput(session, "frameryba", value = input$frameryba + 10)
})

### HYBATKO start
observeEvent(input$plot_click, {
  if (is.null(input$plot_click)) {
    val <- 1
  } else {
    datadir <- datadirInput()
    dataset <- datadir[datadir$rybicka == input$datafish,]
    val <- pow2(nearPoints(dataset, input$plot_click, xvar="head_x", yvar="head_y", addDist = FALSE, maxpoints=1, threshold=100))  #####  !!!!!!!!1
  }
  updateSliderInput(session, "frameryba", value = val)
})

observeEvent(input$loop_video, {
  if (input$frameryba == "") {
    updateSliderInput(session, "frameryba", value = 1)
  }

  if (input$datafish != "all") {
    if (as.integer(input$datafish) < 10) {
      wellnr <- paste("0", input$datafish, sep="")
    } else {
      wellnr <- input$datafish
    }
    write(paste0(input$frameryba,";",input$datadir,";",wellnr,";",input$loop_video[1],";",input$loop_video[2]),"framenr.csv")
  }
})

observeEvent(input$frameryba, {
  if (input$frameryba == "") {
    updateSliderInput(session, "frameryba", value = 1)
  }

  if (input$datafish != "all") {
    if (as.integer(input$datafish) < 10) {
      wellnr <- paste("0", input$datafish, sep="")
    } else {
      wellnr <- input$datafish
    }
    write(paste0(input$frameryba,";",input$datadir,";",wellnr,";",input$loop_video[1],";",input$loop_video[2]),"framenr.csv")
  }

  if ( ! file.exists("browse_running")) {
    if (input$datafish != "all" && input$frameryba != 1) {
      # on windows python
      command <- paste0("python browse_inferred_tar.py ", input$datadir, "/" , wellnr , "/", wellnr, "_*_inferred.csv ", input$datadir, "/" , wellnr , "/images.tar framenr.csv")
      # on linux python3
      #command <- paste0("python3 browse_inferred_tar.py ", input$datadir, "/" , wellnr , "/", wellnr, "_*_inferred.csv ", input$datadir, "/" , wellnr , "/images.tar framenr.csv")
      print(command)
      system(command, wait=F, show.output.on.console=T, invisible=F)
    }
  }

#################################
## DISTPLOT - HEAD_X vs HEAD_Y ##

  output$distPlot <- renderPlot({
  datadir <- datadirInput()

  dataset <- datadir[datadir$rybicka == input$datafish,]

  colors <- c("lightsteelblue3", "moccasin", "lightpink", "lightpink3", "indianred1", "brown1", "red", "red3" , "darkred", "black")
  cuts <- c(0,5,7,10,15,20,35,50,60,70,100)
  bbb <- cuts[1:10]

  hx <- dataset$head_x   # hodnota x
  hy <- dataset$head_y   # hodnota y

  par(mar=c(0,0,0,0), oma=c(0,0,0,0))
  plot(c(0,pixels), c(0,pixels), xlim=c(0,pixels), ylim=c(0,pixels), col="white", axes = FALSE,ann=FALSE)
  par(new=T)

  abline(v=bins_list, col="lightgray", lwd=1.5)
  abline(h=bins_list, col="lightgray", lwd=1.5)
  par(new=T)

  for (i in 1:10) {

    hx1 <- dataset$head_x[dataset$perc <= cuts[i+1] & dataset$perc > cuts[i]]
    hy1 <- dataset$head_y[dataset$perc <= cuts[i+1] & dataset$perc > cuts[i]]

    plot(hx1, hy1, main="", xlab="", ylab="", xlim=c(0,pixels), ylim=c(0,pixels), col=colors[i], pch=19, cex=1.5)
    par(new=T)

  }

  abline(v=0, col="blue4", lwd=5)
  abline(v=pixels, col="blue4", lwd=5)
  abline(h=0, col="blue4", lwd=5)
  abline(h=pixels, col="blue4", lwd=5)

  abline(v=pow3(input$frameryba,dataset), col="chartreuse", lwd=2.5)
  abline(h=pow4(input$frameryba,dataset), col="chartreuse", lwd=2.5)

  par(xpd=TRUE)
  legend(220, 240, bbb, cex=1.0, col=colors, pch=21, lty=1:3)
  })

})
#########

### HYBATKO end

observeEvent(input$lowx, {
  updateSliderInput(session, "delrange", value = c(input$lowx+1,input$highx+1),  min = input$lowx+1, max = input$highx+1, step = 1)
})

observeEvent(input$highx, {
  updateSliderInput(session, "delrange", value = c(input$lowx+1,input$highx+1),  min = input$lowx+1, max = input$highx+1, step = 1)
})

observeEvent({event_data("plotly_relayout", source = "A")}, {
  runjs("var plot = document.getElementById('veloPlot'); Shiny.setInputValue('lowx', parseInt(plot.layout.xaxis.range[0]), {priority: 'event'}); Shiny.setInputValue('highx', parseInt(plot.layout.xaxis.range[1]), {priority: 'event'}); Shiny.setInputValue('lowy', parseInt(plot.layout.yaxis.range[0]), {priority: 'event'}); Shiny.setInputValue('highy', parseInt(plot.layout.yaxis.range[1]), {priority: 'event'});")

   plotlyProxy("veloPlot2",session) %>%
    plotlyProxyInvoke("relayout", list(shapes = list(
         list(type = "rect",
         fillcolor = "blue", line = list(color = "orange"), opacity = 0.3,   # blue
         x0 = input$delrange[1], x1 = input$delrange[2], xref = "x",
         y0 = -5, y1 = input$highy, yref = "y"))))
})

observeEvent({event_data("plotly_click", source = "A")}, {
  event_data_values <- event_data("plotly_click", source = "A")
  updateSliderInput(session, "frameryba", value = event_data_values$x)
})

##############
## VELOPLOT ##
output$veloPlot <- renderPlotly({

  datadir <- datadirInput()
  dataset <- datadir[datadir$rybicka == input$datafish,]
  dataset$velocities[dataset$velocities < input$cutoff_noise] <- input$cutoff_noise

  if (nrow(dataset) != 0) {
    fig <- plot_ly(dataset, x=~X, y=~lengths, type="scatter", mode="lines", line = list(color = "red"), source="A", name="h-t dist")
    fig <- fig %>% add_trace(y=~velocities, mode="lines", line = list(color = "blue"), name="velocity")
    fig <- fig %>% add_trace(y=~acceleration, mode="lines", line = list(color = "magenta"), name="acceleration")
    fig <- fig %>% layout(showlegend = TRUE, yaxis = list(hoverformat = '.1f'), xaxis = list(hoverformat = 'd')) %>%
    event_register('plotly_relayout') %>%
    event_register('plotly_click')
  } else {
    faked <- rep(1,lines_file)
    fig <- plot_ly(dataset, x=1:lines_file, y=faked, type="scatter", mode="lines", line = list(color = "gray"), name="h-t dist del") %>%
    layout(xaxis = list(range=c(1,lines_file))) %>%
    layout(yaxis = list(range=c(-10,50)))
  }
})

########################
## VELOPLOT UNCHANGED ##
output$veloPlot2 <- renderPlotly({

  datadir <- datadirInput()
  dataset <- datadir[datadir$rybicka == input$datafish,]
  lines_file <- nrow(dataset)

  fig <- plot_ly(dataset, x=~X, y=~lengths, type="scatter", mode="lines", line = list(color = "red"))
  fig <- fig %>% add_trace(y=~velocities, mode="lines", line = list(color = "blue"))
  fig <- fig %>% layout(showlegend = FALSE,  yaxis = list(hoverformat = '.1f'), xaxis = list(hoverformat = 'd')) %>%
  add_segments(x = input$frameryba, xend = input$frameryba, y = 0, yend = 50, line = list(color = "chartreuse", width = 2.5)) %>%
  layout(xaxis = list(range=c(1,lines_file))) %>%
  layout(yaxis = list(range=c(-5,50)))
})

#####################################
## TABLES AFTER CHANGING OF A FISH ##
observeEvent(input$datafish, {

  all_data <- datadirInput()
  datadir_after_cutoff <- all_data[all_data$velocities > input$cutoff_noise,]  ## je tu uz cutoff??

  if (input$datafish != "all") {

    rv$data <- as.data.frame(MEAN_LEN_3(input$datadir,input$datafish,datadirInput(),input$cutoff_noise))     ## DATADIR AFTER CUTOFF !!
    rv$deletedRows <- NULL
    rv$deletedRowIndices = list()

    output$MEAN_LEN <- DT::renderDataTable(MEAN_LEN_2(input$datadir,input$datafish,all_data,input$cutoff_noise),selection = 'single')  ## all_data  ????
    output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,all_data,input$cutoff_noise),selection = 'single')  ## PRECO DATADIR_AFTER_CUTOFF ????
    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,all_data,input$cutoff_noise),selection = 'single')  ## all_data   ????
    output$SPEEDUPS_SLOWDOWNS <- DT::renderDataTable(get_speedups_and_slowdowns(input$datadir,input$datafish,all_data,input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise),selection = 'single')  ## ALL DATA


   ## table with DELETES
#    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv"))) {
#      data <- del_table2(input$datadir,input$datafish)
#      if (nrow(data) != 0) {
#        output$MEAN_LEN3 <- DT::renderDataTable(
#       # Add the delete button column
#        deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)
#        )
#      } else {
#        output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
#      }
#    } else {
#      output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
#    }
   ## table with DELETES


    fish_data_len <- nrow(all_data[all_data$rybicka == input$datafish,])
    updateSliderInput(session, "frameryba", max = fish_data_len, value = 1)
    updateSliderInput(session, "delrange", max = fish_data_len, value = c(1,fish_data_len))


    ## update cutoff according to a fish  (each fish can have its own set of cutoffs)
    table_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)
    table_for_fish <- table_csv[table_csv$directory == input$datadir & table_csv$fish == input$datafish,]

    updateNumericInput(session, "cutoff_speedup", "cutoff_speedup", value = table_for_fish$cutoff_speedup)
    updateNumericInput(session, "skip_speedup", "skip_speedup", value = table_for_fish$skip_speedup)
    updateNumericInput(session, "cutoff_slowdown", "cutoff_slowdown", value = table_for_fish$cutoff_slowdown)
    updateNumericInput(session, "skip_slowdown", "skip_slowdown", value = table_for_fish$skip_slowdown)
    updateNumericInput(session, "cutoff_noise", "cutoff_velocity_noise", value = table_for_fish$cutoff_noise)

    output$info2 <- renderPrint({
      data_mov <- get_number_of_movement_frames(input$datadir,input$datafish,datadirInput(),input$cutoff_speedup,input$cutoff_count_movements,input$cutoff_noise)
      print(paste0(data_mov[1], " frames (", data_mov[2], " %) moving"))
    })

    output$info3 <- renderPrint({
      perc_matrix <- get_mean_percents_selected_areas(input$datadir,input$datafish, datadirInput())
      areas_all <- get_mean_per_area(perc_matrix)

      area_last <- tail(areas_all,1)
      for (area_num in 0:area_last) {
        pos <- area_num*2+1
        mean_per_area <- areas_all[pos]

        print(paste0(round(mean_per_area,digits=1), " % in area ", area_num))
      }
    })

    output$info4 <- renderPrint({
      maskalist <- get_maskalist_selected_areas(input$datadir,input$datafish, datadirInput())
      crossings_list <- strsplit(get_crossings(maskalist), ';')
      len_cross <- length(crossings_list)
      for (num_cross in 1:len_cross) {
#        if (num_cross == 1) {
#          crossings_output <- crossings_list[num_cross]
#        } else {
#          crossings_output <- rbind(crossings_output,crossings_list[num_cross])
#        }
#        print(crossings_list[[1]][num_cross])
        print(crossings_list[[num_cross]])
      }
#      print(crossings_output)
    })

    output$maskaplot <- renderPlot({
      selected_areas <- as.matrix(read.csv("ryby_maska.csv", sep=",", header=FALSE))
      par(mar=c(0,0,0,0), oma=c(0,0,0,0))
      image(t(apply(selected_areas, 2, rev)), axes=F)
      #heatmap(selected_areas, Rowv = NA, Colv = "Rowv")
    })

  } else {

    # set default values from header of this file
    updateNumericInput(session, "cutoff_speedup", "cutoff_speedup", value = cutoff_speedup)
    updateNumericInput(session, "skip_speedup", "skip_speedup", value = skip_speedup)
    updateNumericInput(session, "cutoff_slowdown", "cutoff_slowdown", value = cutoff_slowdown)
    updateNumericInput(session, "skip_slowdown", "skip_slowdown", value = skip_slowdown)
    updateNumericInput(session, "cutoff_noise", "cutoff_velocity_noise", value = cutoff_noise)

  }

})

observeEvent(input$datadir, {

  all_data <- datadirInput()
  datadir_after_cutoff <- all_data[all_data$velocities > input$cutoff_noise,]  ## je tu uz cutoff??

  if (input$datafish != "all") {

    rv$data <- as.data.frame(MEAN_LEN_3(input$datadir,input$datafish,datadirInput(),input$cutoff_noise))     ## DATADIR AFTER CUTOFF !!
    rv$deletedRows <- NULL
    rv$deletedRowIndices = list()

    output$MEAN_LEN <- DT::renderDataTable(MEAN_LEN_2(input$datadir,input$datafish,all_data,input$cutoff_noise),selection = 'single')  ## all_data  ????
    output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,all_data,input$cutoff_noise),selection = 'single')  ## PRECO DATADIR_AFTER_CUTOFF ????
    output$MEAN_LEN4 <- DT::renderDataTable(MEAN_LEN_4(input$datadir,input$datafish,all_data,input$cutoff_noise),selection = 'single')  ## all_data   ????
    output$SPEEDUPS_SLOWDOWNS <- DT::renderDataTable(get_speedups_and_slowdowns(input$datadir,input$datafish,all_data,input$cutoff_speedup,input$skip_speedup,input$skip_slowdown,input$cutoff_noise),selection = 'single')  ## ALL DATA


   ## table with DELETES
#    if (file.exists(paste0("deletes_", input$datadir, "_", input$datafish,".tsv"))) {
#      data <- del_table2(input$datadir,input$datafish)
#      if (nrow(data) != 0) {
#        output$MEAN_LEN3 <- DT::renderDataTable(
#       # Add the delete button column
#        deleteButtonColumn(rv$data, 'delete_button',input$datafish,input$datadir,input$cutoff_noise)
#        )
#      } else {
#        output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
#      }
#    } else {
#      output$MEAN_LEN3 <- DT::renderDataTable(MEAN_LEN_3(input$datadir,input$datafish,datadir_orig,input$cutoff_noise),selection = 'single')
#    }
   ## table with DELETES


    fish_data_len <- nrow(all_data[all_data$rybicka == input$datafish,])
    updateSliderInput(session, "frameryba", max = fish_data_len, value = 1)
    updateSliderInput(session, "delrange", max = fish_data_len, value = c(1,fish_data_len))

    ## update cutoff according to a fish  (each fish can have its own set of cutoffs)
    table_csv <- read.csv("RESULTS_conditions.csv", sep=",", header=TRUE, row.names = 1)
    table_for_fish <- table_csv[table_csv$directory == input$datadir & table_csv$fish == input$datafish,]

    updateNumericInput(session, "cutoff_speedup", "cutoff_speedup", value = table_for_fish$cutoff_speedup)
    updateNumericInput(session, "skip_speedup", "skip_speedup", value = table_for_fish$skip_speedup)
    updateNumericInput(session, "cutoff_slowdown", "cutoff_slowdown", value = table_for_fish$cutoff_slowdown)
    updateNumericInput(session, "skip_slowdown", "skip_slowdown", value = table_for_fish$skip_slowdown)
    updateNumericInput(session, "cutoff_noise", "cutoff_velocity_noise", value = table_for_fish$cutoff_noise)

    output$info2 <- renderPrint({
      data_mov <- get_number_of_movement_frames(input$datadir,input$datafish,datadirInput(),input$cutoff_speedup,input$cutoff_count_movements,input$cutoff_noise)
      print(paste0(data_mov[1], " frames (", data_mov[2], " %) moving"))
    })


    output$info3 <- renderPrint({
      perc_matrix <- get_mean_percents_selected_areas(input$datadir,input$datafish, datadirInput())
      areas_all <- get_mean_per_area(perc_matrix)

      area_last <- tail(areas_all,1)
      for (area_num in 0:area_last) {
        pos <- area_num*2+1
        mean_per_area <- areas_all[pos]

        print(paste0(round(mean_per_area,digits=1), " % in area ", area_num))
      }
    })

    output$info4 <- renderPrint({
      maskalist <- get_maskalist_selected_areas(input$datadir,input$datafish, datadirInput())
      crossings_list <- strsplit(get_crossings(maskalist), ';')
      len_cross <- length(crossings_list)
      for (num_cross in 1:len_cross) {
#        if (num_cross == 1) {
#          crossings_output <- crossings_list[num_cross]
#        } else {
#          crossings_output <- rbind(crossings_output,crossings_list[num_cross])
#        }
#        print(crossings_list[[1]][num_cross])
        print(crossings_list[[num_cross]])
      }
#      print(crossings_output)
    })

  } else {

    # set default values from header of this file
    updateNumericInput(session, "cutoff_speedup", "cutoff_speedup", value = cutoff_speedup)
    updateNumericInput(session, "skip_speedup", "skip_speedup", value = skip_speedup)
    updateNumericInput(session, "cutoff_slowdown", "cutoff_slowdown", value = cutoff_slowdown)
    updateNumericInput(session, "skip_slowdown", "skip_slowdown", value = skip_slowdown)
    updateNumericInput(session, "cutoff_noise", "cutoff_velocity_noise", value = cutoff_noise)
  }
})

##########################
## vypisat text do boxu ##
  output$info <- renderPrint({

  datadir2 <- datadirInput()
  dataset <- datadir2[datadir2$rybicka == input$datafish,]

  if (as.integer(input$datafish) < 10) {
  fish <- paste("0", input$datafish, sep="")
  } else {
  fish <- input$datafish
  }
    show_vars <- c(input$show_vars1, input$show_vars2, input$show_vars3, input$show_vars4, input$show_vars5)
    print(show_vars)

  })

  output$info2 <- renderPrint({
    print(input$cutoff_velo)
  })


#############
## SUMMARY ##
  output$summary <- renderPrint({
    datadir <- datadirInput()
    summary(datadir)

  })
}

shinyApp(ui = ui, server = server)
#runApp(list(ui = ui, server = server), launch.browser = TRUE)
