
# Finds the coordinates of the Moore neighbors on a non-toroidal grid.
# Excludes ego.
mooreNeigh <- \(x, y, minx = 1, miny = 1, maxx, maxy) {
  nb <- expand.grid(
    x = max(x - 1, minx):min(x + 1,maxx),
    y = max(y - 1, miny):min(y + 1,maxy)
  )
  
  #remove ego and return:
  return(nb[!(nb$x == x & nb$y == y),])
}
#mooreNeigh(x = 4, y = 1, minx = 1, miny = 1, maxx = 10, maxy = 10) |> plot()


nbs <- \(m, x, y) {
  
  # find moore neighborhood by subsetting matrix around ego
  nb <- m[max(1, y-1) : min(nrow(m), y+1) , max(1, x-1) : min(ncol(m), x+1)]
  
  # exclude ego
  ifelse(
    m[y,x], # If ego is inhabited...
    return(sum(nb) - 1), #... we remove ego from the tally ...
    return(sum(nb)) # else, we return the whole tally of neighbors
  )
}

# Weight calculation for the RI-model with one opinion dimension and one group
# characteristic. Opinions are in [0,1] and groups are {0,1}.
calcW <- \(oi, oj, gi, gj, H = 0.6, mu = 0.5) {
  mu * (1 - (2 * ((abs(oj - oi) * H) + (abs(gj - gi) * (1 - H)))))
}

# Generalized version of weight calculation. 
# mu scales the strength of influence; th sets the opinion distance that
# triggers repulsion.
calcW_gen <- \(oi, oj, mu = 0.5, th = 0.4) { # th must be in (0,1] -> zero excluded
  diff <- abs(oj - oi)
  w <- sapply(diff, FUN = \(diff) ifelse(
    diff <= th,
    w <- mu * (1 - (diff / th)), # assimilation
    w <- -1 * mu * ((diff - th) / (1 - th)) # repulsion
  ))
  return(w)
}


# Rule for one-to-one influence
newO <- \(oi, oj, w) oi + (w * (oj - oi))


# Rule for many-to-one influence.
newO_mto <- \(oi, oNeigh, mu = 0.5, mode = "mean") {
  w <- calcW(oi = oi, oj = oNeigh, mu = mu)
  
  if (mode == "mean") { # interact with an average of the neighbors
    o <- oi + (mean(w, na.rm = TRUE) * (mean(oNeigh, na.rm = TRUE) - oi))
  }
  
  if (mode == "min") { # interact with the most influential neighbor
    alter <- which(abs(w) == min(abs(w), na.rm = TRUE))[1]
    o <- oi + (w[alter] * (oNeigh[alter] - oi))
  }
  
  if (mode == "max") { # interact with the most influential neighbor
    alter <- which(abs(w) == max(abs(w), na.rm = TRUE))[1]
    o <- oi + (w[alter] * (oNeigh[alter] - oi))
  }
  return(truncate(o, min = 0, max = 1))
}


rleDecoder <- \(pat) {
  # separating input into lines
  pat <- strsplit(x = pat, split = "\n")
  
  # concatenating data lines:
  pat <- paste0(pat[[1]][2:length(pat[[1]])], collapse = "")
  
  # removing end character "!":
  pat <- gsub(pattern = "!", replacement = "", x = pat)
  
  # Finding which instances we need to repeat -- we do this by identifying what
  # integer numbers there are in the pat string:
  digits <- unlist(regmatches(pat, gregexpr("[[:digit:]]+", pat)))
  
  # Essentially, for each integer d we look up where it occurs first, and we
  # repeat d times the character that follows it.
  for (d in 1:length(digits)) {
    
    # where the d occurs first:
    i = stringr::str_locate(pat, digits[d])[1,]
    
    # what character is to be repeated:
    char <-  substr(pat, start = i[2] + 1, stop =  i[2] + 1)
    
    # this is the expanded string with the repeated character:
    expanded <- paste0(rep(char, times = as.numeric(digits[d])), collapse = "")
    
    # erasing the part of the string that we are replacing:
    #stringr::str_sub(pat, start = i[1], end = i[2] + 1) <- ""
    
    # injecting the expanded string at the needed location:
    ifelse(
      i[1] == 1,
      patStart <- "",
      patStart <- stringr::str_sub(pat, start = 1, end = i[1] - 1)
    )
    ifelse(
      i[2] == nchar(pat),
      patEnd <- "",
      patEnd <- stringr::str_sub(pat, start = i[2] + 2, end = nchar(pat))
    )
    pat <- paste0(patStart, expanded, patEnd, collapse = "")
  }
  
  # character "$" signifies end-of-line. Therefore:
  pat <- strsplit(x = pat, split = "\\$")[[1]]
  
  # creating a matrix of the necessary size:
  patm <- matrix(0, nrow = length(pat), ncol = max(nchar(pat)))
  
  # filling in the matrix, row by row:
  for (i in 1:length(pat)) {
    ones <- stringr::str_locate_all(string = pat[i], pattern = "o")[[1]][,1]
    patm[i, ones] <- 1
  }
  
  return(patm)
}
#rleDecoder(pat = "x = 34, y = 31, rule = B3/S23
#33bo$31b2o$32b2o$9bo$bo8bo$2bo5b3o$3o3$5bo$6bo$4b3o$24bobo$25b2o$25bo
#2$27bobo$27b2o$28bo$31b3o$31bo$32bo7$5b2o$6b2o$5bo!") |> image()


calcDist <- function(x1, x2, y1, y2, fastMode = TRUE) {
  ifelse(
    fastMode,
    return(((x2 - x1)^2) + ((y2 - y1)^2)),
    return(sqrt(((x2 - x1)^2) + ((y2 - y1)^2)))
  )
}


################################################################################
# General purpose functions

normalize <- \(x) (x - min(x)) / (max(x) - min(x))

truncate <- \(x, min = 0, max = 1) {
  x[x < min] <- min
  x[x > max] <- max
  return(x)
}