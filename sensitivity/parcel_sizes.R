mean_parcel <- log(10)
sd_parcel <- log(2)

n_acres <- ncell*2224

# generate parcel sizes from lognormal distribution
parcel_sizes <- rlnorm(n_acres/exp(mean_parcel), mean_parcel, sd_parcel)
parcel_sizes <- parcel_sizes[cumsum(parcel_sizes)<=n_acres]

# assign parcels to grid cells
id.df <- data.frame(grid=1, parcel=1, parcel.acres=parcel_sizes[1])
j <- 2 # parcel index
r <- 2 # row index
for(i in 1:10) { # grid cell index
  while(sum(id.df$parcel.acres[id.df$grid==i]) < 20) {
    id.df <- rbind(id.df, 
                   data.frame(grid=i, parcel=j, parcel.acres=parcel_sizes[j]))
    r <- r+1
    j <- j+1
  }
  while(sum(id.df$parcel.acres[id.df$grid==i]) > 20) {
    over <- sum(id.df$parcel.acres[id.df$grid==i]) - 20
    if(over > 20) {
      over_orig <- over
      over <- over_orig - 20
    }
    ## partition problem??
    
    while(over > 0) {
      if(over > 20) {
        overover <- over - 20
        over <- over - 20
      }
      id.df$parcel.acres[r-1] <- id.df$parcel.acres[r-1] - over
      id.df <- rbind(id.df, 
                     data.frame(grid=i+1, 
                                parcel=j-1, 
                                parcel.acres=over))
      i <- i+1
      r <- r + 1
    }
  }
}

