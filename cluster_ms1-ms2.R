###############################################################################
# AUTH  Jeff Jones
# DATE  2019.01.19
#
# 

options(warn=-1)

###############################################################################
# USER INPUT
outf <- x_df <- y_df          <- NULL
x_mz_col <- y_mz_col          <- 'mz'
x_rt_col <- y_rt_col          <- 'rt'
x_id_col <- y_id_col          <- 'id'
tol_mz_da                     <- 0.1
tol_rt_sec                    <- 10 
segment_size                  <- 1000

for (arg in commandArgs()){
  arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
  if( grepl("--out", arg) ) outf <- arg_value
  if( grepl("--nf", arg) ) x_df <- arg_value
  if( grepl("--mf", arg) ) y_df <- arg_value
  if( grepl("--xmz", arg) ) x_mz_col <- arg_value
  if( grepl("--xrt", arg) ) x_rt_col <- arg_value
  if( grepl("--xid", arg) ) x_id_col <- arg_value
  if( grepl("--ymz", arg) ) y_mz_col <- arg_value
  if( grepl("--yrt", arg) ) y_rt_col <- arg_value
  if( grepl("--yid", arg) ) y_id_col <- arg_value
  if( grepl("--tmz", arg) ) tol_mz_da <- arg_value
  if( grepl("--trt", arg) ) tol_rt_sec <- arg_value
  if( grepl("--seg", arg) ) segment_size <- arg_value
}

###############################################################################
# INPUT VALIDATION
message <- NULL
if(is.null(outf)) message <- paste0(message, "  no output file declared\n")
if(is.null(x_df)) message <- paste0(message, "  no reference file declared\n")
if(is.null(y_df)) message <- paste0(message, "  no cluster file declared\n")
if(!file.exists(x_df)) message <- paste0(message, "  reference file (--nf) not found\n")
if(!file.exists(y_df)) message <- paste0(message, "  cluster file (--mf) not found\n")
if(!grepl("rds|csv$", x_df)) message <- paste0(message, "  reference file (--mf) not a supported format\n")
if(!grepl("rds|csv$", y_df)) message <- paste0(message, "  cluster file (--mf) not a supported format\n")
if(!grepl("rds|csv$", outf)) message <- paste0(message, "  output file (--out) not a supported format\n")

if(!is.null(message)) stop("ERROR\n", message)

cat("cluster ms1-ms2 started\n")
cat(" reference file:      ", x_df, "\n")
cat(" cluster file:        ", y_df, "\n")
cat(" output file:         ", outf, "\n")


###############################################################################
# FUNCTION: input validation on the file type
read.input <- function(file_name){
  if(grepl(".rds$", file_name)){
    df <- readRDS(file_name)
  } else if(grepl(".csv$", file_name)){
    df <- read.csv(file_name)
  }
  return(df)
}

###############################################################################
# FUNCTION: calculate all possible differences
apdif <- function(x_df=c(),
                  y_df=c(),
                  x_col_val=NULL,
                  x_col_ids=NULL,
                  y_col_val=NULL,
                  y_col_ids=NULL){
  
  if(!x_col_val %in% colnames(x_df))
    stop("\n\t=> '", x_col_val , "' not a valid column in data.frame")
  
  if(!y_col_val %in% colnames(y_df))
    stop("\n\t=> '", y_col_val , "' not a valid column in data.frame")
  
  if(!x_col_ids %in% colnames(x_df))
    stop("\n\t=> '", x_col_ids , "' not a valid column in data.frame")
  
  if(!y_col_ids %in% colnames(y_df))
    stop("\n\t=> '", y_col_ids , "' not a valid column in data.frame")
  
  
  
  
  x_val                     <- x_df[,x_col_val]
  y_val                     <- y_df[,y_col_val]
  
  x_ids                      <- x_df[,x_col_ids]
  y_ids                      <- y_df[,y_col_ids]
  
  x_n                        <- length(x_val)
  y_n                        <- length(y_val)
  
  x_m                        <- array(x_val, c(x_n, y_n))
  
  all_dif                    <- as.numeric(t(x_m) - y_val)
  
  x_m_ids                    <- as.character(t(array(x_ids, c(x_n, y_n))))
  y_m_ids                    <- as.character(array(y_ids, c(y_n, x_n)))
  
  df <- data.frame(dif=all_dif,
                   row_x=x_m_ids,
                   row_y=y_m_ids)
  
  # rename the columns
  
  names(df)[2] <- x_col_ids
  names(df)[3] <- y_col_ids
  
  if(x_col_ids == y_col_ids){
    names(df)[2] <- paste0(names(df)[2], "_x")
    names(df)[3] <- paste0(names(df)[3], "_y")
  }
  
  return(df)
  
}

###############################################################################
# FUNCTION: cluster based on a 2d value coordinate
library(progress)
cluster <- function(x_df,
                    y_df,
                    x_mz_col = 'mz',
                    y_mz_col = 'mz',
                    x_rt_col = 'rt',
                    y_rt_col = 'rt',
                    x_id_col = 'id',
                    y_id_col = 'id',
                    tol_mz_da = 0.1,
                    tol_rt_sec = 10,
                    segment_size = 2000){
  
  cluster_df <- c()
  
  # determine the N number of segments
  n_ydf <- dim(y_df)[1]
  n_seg <- ceiling(n_ydf/segment_size)
  
  x_df <- x_df[order(x_df[,x_mz_col]),]
  y_df <- y_df[order(y_df[,y_mz_col]),]
  
  pb <- progress_bar$new(
    format = " clustering [:bar] :percent eta: :eta",
    total = n_seg, clear = FALSE, width= 60)
  
  for( i in 1:n_seg ){
    
    pb$tick()
    
    i_min <- max(1, (i-1) * segment_size)
    i_max <- min(i * segment_size, n_ydf)
    
    this_y_df <- y_df[i_min:i_max,]
    
    x_mz_min <- min(this_y_df[,y_mz_col]) - tol_mz_da * 2
    x_mz_max <- max(this_y_df[,y_mz_col]) + tol_mz_da * 2
    x_rt_min <- min(this_y_df[,y_rt_col]) - tol_rt_sec * 2
    x_rt_max <- max(this_y_df[,y_rt_col]) + tol_rt_sec * 2
    
    this_x_df <- x_df
    this_x_df <- this_x_df[this_x_df[,x_mz_col] >= x_mz_min,]
    this_x_df <- this_x_df[this_x_df[,x_mz_col] < x_mz_max,]
    this_x_df <- this_x_df[this_x_df[,x_rt_col] >= x_rt_min,]
    this_x_df <- this_x_df[this_x_df[,x_rt_col] < x_rt_max,]
    
    if(dim(this_x_df)[1] == 0)
      next()
    
    delta_mz_df <- apdif(this_x_df, this_y_df,
                         x_mz_col, x_id_col,
                         y_mz_col, y_id_col)
    
    delta_lc_df <- apdif(this_x_df, this_y_df,
                         x_rt_col, x_id_col,
                         y_rt_col, y_id_col)
    
    delta_mz_df <- delta_mz_df[abs(delta_mz_df$dif) <= tol_mz_da,]
    delta_lc_df <- delta_lc_df[abs(delta_lc_df$dif) <= tol_rt_sec,]
    
    m_x <- x_id_col
    m_y <- y_id_col
    
    if(x_id_col == y_id_col){
      m_x <- paste0(m_x, "_x")
      m_y <- paste0(m_y, "_y")
    }
    
    this_cluster_df <- merge(delta_mz_df, delta_lc_df, by=c(m_x, m_y), suffixes=c("_mz", "_lc"))
    
    cluster_df <- rbind(cluster_df, this_cluster_df)
  }
  
  return(cluster_df)
  
}


x_df <- read.input(x_df)
y_df <- read.input(y_df)

out_df <- cluster(x_df,
               y_df,
               x_mz_col,
               y_mz_col,
               x_rt_col,
               y_rt_col,
               x_id_col,
               y_id_col,
               tol_mz_da,
               tol_rt_sec,
               segment_size
)


if(grepl("csv$", outf))
  write.csv(out_df, outf)

if(grepl("rds$", outf))
  saveRDS(out_df, outf)

cat("\nFINISHED: data written to ", outf, "\n")

