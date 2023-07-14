# Latest update: Thomas Feliciani, 14 July 2023 
# This script reproduces the artwork images for a manuscript in preparation.
# Runs in R 4.3.1



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Loading resources

rm(list = ls())

library("compiler")
library("reshape2")
library("dplyr")
library("lhs")
library("png")
library("ambient")
#library("viridis")
library("plotly")
library("igraph")
library("ggplot2")
library("rayshader")
library("rgl")

source("util.r")





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Cover image

seed = 12345
set.seed(seed)


createBlock <- \(
  lots = 100,
  density = 0.8,
  probMatchHeight = 0.3,
  maxHeight = 1
) {
  m <- matrix(rep(0, times = lots), nrow = round(sqrt(lots)))
  
  for (i in 1:(lots * density)) {
    vacant <- TRUE
    while (vacant) {
      # Pick a random grid cell.
      x <- sample(1:round(sqrt(lots)), size = 1)
      y <- sample(1:round(sqrt(lots)), size = 1)
      
      #print(i)
      # if we found an empty grid cell...
      if (m[y,x] == 0) {
        vacant <- FALSE
        
        # Then place check for adjacent pre-existing buildings.
        nbHeights <- c(
          m[min(round(sqrt(lots)), y + 1), x], # N
          m[max(1, y - 1), x], # S
          m[y, min(round(sqrt(lots)), x + 1)], # E
          m[y, max(1, x - 1)] # W
        )
        
        # If there are no building nearby create a new randomly-tall one. 
        # Else, create a new one that matches the shortest neighbor:
        ifelse(
          any(nbHeights != 0) & runif(1) < probMatchHeight,
          m[y, x] <- min(nbHeights[nbHeights > 0]),
          #m[y, x] <- runif(n = 1, min = 0.1, max = maxHeight)
          m[y, x] <- sample(seq(from = 0.1, to = maxHeight, by = 0.1), size = 1)
        )
      }
    }
  }
  
  block <- reshape2::melt(m)
  names(block) <- c("y", "x", "z")
  return(block[block$z != 0,])
}
#block <- createBlock()


mGroup <- matrix(0, nrow = 6, ncol = 6)
mGroup[1, 2:6] <- 1
mGroup[2, 3:6] <- 1
mGroup[3, 3:6] <- 1
mGroup[4, 4:6] <- 1
mGroup[5, 4:6] <- 1
mGroup[6, 5:6] <- 1

mPol <- matrix(5, nrow = 6, ncol = 6)
mPol[1, 1:2] <- mPol[2, 2:3] <- mPol[3, 2:3] <- mPol[4, 3:4] <-
  mPol[5, 3:4] <- mPol[6, 4:5] <- 1 # Boundary blocks
mPol[1, 3] <- mPol[2, 1] <- mPol[2, 4] <- mPol[3, 1] <- mPol[3, 4] <-
  mPol[4, 2] <- mPol[4, 5] <- mPol[5, 2] <- mPol[5, 5] <- mPol[6, 3] <- 
  mPol[6, 6] <- 2 # Distance 2 from border
mPol[1,4] <- mPol[2,5] <- mPol[3,5] <- mPol[4,1] <- mPol[4,6] <- mPol[5,1] <-
  mPol[5,6] <- mPol[6,2] <- 3 # Distance 3
mPol[1,5] <- mPol[2,6] <- mPol[3,6] <- mPol[6,1] <- 4 # Distance 4
dPol <- reshape2::melt(mPol)
names(dPol) <- c("r", "c", "col")

mDens <- matrix(0.25, nrow = 6, ncol = 6)
mDens[2:5, 2:5] <- 0.4
mDens[3:4, 3:4] <- 0.6 # downtown
mDens[4,3] <- 0.7 # center
mDens[1,1] <- mDens[1,6] <- mDens[6,1] <- mDens[6,6] <- 0.15 # corners
mDens[2,2] <- mDens[2,5] <- mDens[5,2] <- mDens[5,5] <- 0.32
#mDens[4,4] <- 0.55

for (r in 1:6) for (c in 1:6) {
  #print(paste(r, c))
  block <- createBlock(
    lots = 100,
    density = mDens[r,c],
    probMatchHeight = 0.3,
    maxHeight = (mDens[r,c] + 0.1) * 3
  )
  block$group <- mGroup[r,c]
  block$r <- r
  block$c <- c
  
  ifelse(
    r == 1 & c == 1,
    d <- block,
    d <- rbind(d, block)
  )
}

ggtheme <- ggplot2::theme(
  strip.background = element_blank(),
  strip.text = element_blank(),
  legend.position = "NA",
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  panel.spacing = unit(3, "pt")
)

tileCols <- gray.colors(n = max(dPol$col) + 2, start = 0.65, end = 0.95)
tileCols <- tileCols[-c(2,3)]

ggp <- ggplot(
  data = d, 
) +
  geom_rect( # Facet background at distance 1
    data = dPol[dPol$col == 1,], fill = tileCols[1],
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_rect( # Facet background at distance 2
    data = dPol[dPol$col == 2,], fill = tileCols[2],
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_rect( # Facet background at distance 3
    data = dPol[dPol$col == 3,], fill = tileCols[3],
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_rect( # Facet background at distance 4
    data = dPol[dPol$col == 4,], fill = tileCols[4],
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_rect( # Facet background at distance 5
    data = dPol[dPol$col == 5,], fill = tileCols[5],
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  geom_tile(
    aes(x = x, y = y, fill = as.factor(group), color = as.factor(group))#,
    #linewidth = 1
  ) +
  facet_grid(r ~c) + 
  scale_fill_manual(values = c("1" = "orange", "0" = "#EEEEEE")) +
  scale_color_manual(values = c("1" = "orange", "0" = "#EEEEEE")) +
  scale_x_continuous(limits = c(1, 10), expand = c(0.06,0.06)) +
  scale_y_continuous(limits = c(1, 10), expand = c(0.06,0.06)) +
  ggtheme
  
ggph <- ggplot(data = d, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  facet_grid(r ~c) +
  scale_x_continuous(limits = c(1, 10), expand = c(0.06,0.06)) +
  scale_y_continuous(limits = c(1, 10), expand = c(0.06,0.06)) +
  ggtheme

print(ggp)

rayshader::plot_gg(
  ggobj = ggp,
  ggobj_height = ggph,
  background = "white",#"#FFFFFF00",
  height = 4,
  width = 4,
  scale = 50,#30, # scales the vertical extrusion
  offset_edges = FALSE,
  multicore = TRUE, 
  height_aes = "fill",
  sunangle = 210, # 315 = NW; 270 = W
  shadow_intensity = 0.3,
  fov = 70,
  phi = 23, # camera altitude
  theta = 45, # 45 = from SE
  zoom = 0.45 # higher is farther
)
#render_depth(preview_focus = TRUE)
#Focus distance: 1774.77
#Focal range: 1139.33-2792.74

rayshader::render_depth(
  #focus = 1659.06, 
  focallength = 700,
  filename = "./output/cover.png",
  software_render = TRUE, # allows higher resolution than the rgl window
  height = 1800,
  width = 3000,
  camera_location = c(1155.12, 650, 1155.12), #c(1155.12, 693.42, 1155.12),
  camera_lookat = c(0, -100, 0), #c(0, 0, 0),
  print_scene_info = TRUE
)

# clean enviornment:
rgl::close3d()
rm(
  block, d, dPol, ggp, ggph, ggtheme, mDens, mGroup, mPol, 
  c, r, seed, tileCols, createBlock
)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Chapter 1 - Introduction and summary
#
# ______________________________________________________________________________
# Chapter cover:
# Tree

set.seed(5)

#genTree <- function(x, y, a, d, parent) {
genTree <- function(x, y, a, angleNoise = 12, d, col = 0, gen = 1) {
  x2 <- y2 <- 0
  #d2r <- pi / 180 
  a1 <- a * (pi / 180)
  d1 = 0
  if (d <= 0) return()
  if (d > 0) {
    d1 = d * 10
    x2 = x + cos(a1) * d1
    y2 = y + sin(a1) * d1
    
    
    branch <- c(
      x0 = x,
      y0 = y,
      x1 = x2,
      y1 = y2,
      col = col,
      gen = gen
    )
    
    tree[nrow(tree) + 1,] <<- branch
    
    a1 <- a - 20 + rnorm(n = 1, mean = 0, sd = angleNoise)
    a2 <- a + 20 + rnorm(n = 1, mean = 0, sd = angleNoise)
    
    if (col == 0) {
      genTree(x = x2, y = y2, a = a1, d = d - 1, col = -1, gen = gen + 1)
      genTree(x = x2, y = y2, a = a2, d = d - 1, col = 1, gen = gen + 1)
    } else {
      genTree(x = x2, y = y2, a = a1, d = d - 1, col = col * 1.5, gen = gen + 1)
      genTree(x = x2, y = y2, a = a2, d = d - 1, col = col * 1.5, gen = gen + 1)
    }
  }
}

tree <- data.frame(
  x0 = numeric(),
  y0 = numeric(),
  x1 = numeric(),
  y1 = numeric(),
  col = numeric(),
  gen = numeric()
)
genTree(x = 1, y = 1, a = 90, d = 8)

#tree$x0 <- normalize(tree$x0)
#tree$x1 <- normalize(tree$x1)
#tree$y0 <- normalize(tree$y0)
#tree$y1 <- normalize(tree$y1)


shadow <- tree
for (i in 1:nrow(shadow)) {
  
}

shadow$y0 <- - shadow$y0 ** 0.5
shadow$y1 <- - shadow$y1 ** 0.5
shadow$x0 <- shadow$x0 ** 1.1
shadow$x1 <- shadow$x1 ** 1.1

png(
  filename = "./output/Chapter_1.png",
  width = 2200, height = 1600, res = 300, units = "px"
)
ggplot(
  data = tree,
  aes(x = x0, y = y0, xend = x1, yend = y1, color = col, linewidth = gen)
) + 
  geom_segment(data = shadow, color = "gray95") +
  geom_segment() +
  scale_color_gradientn(colors = c(
    rep("orange", times = 3),
    "gray20",
    rep("gray80", times = 3)
  )) +
  scale_linewidth_continuous(range = c(3, 0.5)) +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "NA"
  )
dev.off()

rm(tree, shadow, i, genTree)

# ______________________________________________________________________________
# Figure 1.1. -- attitude polarization and alignment

set.seed(1234567)

d3 <- d2 <- d1 <- data.frame(
  panel = rep(NA, times = 100),
  group = c(rep(1, times = 50), rep(2, times = 50)),
  x = rep(NA, times = 100)
)

# no polarization
d1$panel <- "no polarization" 
d1$x <- rbeta(n = 100, shape1 = 1.5, shape2 = 1.5)

d2$panel <- "polarization"
randsample <- sample(c(TRUE, FALSE), size = 100, replace = TRUE)
d2$x[randsample] <- rbeta(n = sum(randsample), shape1 = 1.5, shape2 = 10)
d2$x[!randsample] <- rbeta(n = sum(!randsample), shape1 = 10, shape2 = 1.5)

d3$panel <- "polarization and alignment"
d3$x[d3$group == 1] <- rbeta(n = 50, shape1 = 1.5, shape2 = 5)
d3$x[d3$group == 2] <- rbeta(n = 50, shape1 = 5, shape2 = 1.5)


png(
  filename = "./output/Fig_1.1.png",
  width = 2000, height = 600, res = 300, units = "px"
)
ggplot(data = rbind(d1, d2, d3), aes(x = x, fill = as.factor(group))) + 
  geom_density(color = NA, alpha = 0.5) +
  facet_wrap(
    ~ as.factor(panel),
    #labeller = as_labeller(c("1" = "attitude", "2" = "attitude", "3" = "attitude"))
  ) +#, scales = "free_y") +
  scale_fill_manual(values = c("gray65", "gray30"), labels = c("", "")) +
  labs(x = "attitude", y = "frequency", fill = "ethnic groups") +
  scale_y_continuous(limits = c(0, 3.2), expand = c(0,0)) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0,0),
    breaks = c(0.1, 0.9),#breaks = c(0, 1),
    labels = c("0" = "-", "1" = "+")
  ) +
  theme(
    panel.background = element_rect(fill = "gray97"),
    plot.background = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    strip.background = element_rect(fill = "gray97"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.position = "left",
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5))
dev.off()

rm(d1, d2, d3, randsample)

# ______________________________________________________________________________
#
# Figure 1.2. -- max polarization or max alignment?
d1 <- data.frame(
  group = c(rep(1, times = 4), rep(2, times = 10)),
  x = c(1, 1.3, 1, 1.3, 1.7, 1.8, 1.9, 2, 2.1, 1.7, 1.8, 1.9, 2, 2.1),
  attitude = c(0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
  y = c(0.2, 0.2, 0.8, 0.8, 0, 0.2, 0, 0.2, 0, 0.8, 1, 0.8, 1, 0.8),
  condition = "max polarization"
)
d2 <- data.frame(
  group = c(rep(1, times = 4), rep(2, times = 10)),
  x = c(1, 1.3, 1, 1.3, 1.7, 1.8, 1.9, 2, 2.1, 1.7, 1.8, 1.9, 2, 2.1),
  attitude = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  y = c(0.2, 0.2, 0.8, 0.8, 0, 0.2, 0, 0.2, 0, 0.8, 1, 0.8, 1, 0.8),
  condition = "max alignment"
)

d <- as.data.frame(rbind(d1, d2))
d$y = d$y * 0.8
#d$x = 1
#d$y <- d$attitude
d$condition <- factor(
  x = d$condition,
  levels = c("max polarization", "max alignment"),
  labels = c("max polarization", "max alignment")
)
d$group <- factor(
  x = d$group,
  levels = c(1, 2),
  labels = c("ethnic group A", "ethnic group B")
)
d$attitude <- factor(
  x = d$attitude,
  levels = c(0, 1),
  labels = c("-", "+")
)

png(
  filename = "./output/Fig_1.2.png",
  width = 1800, height = 900, res = 300, units = "px"
)
ggplot(
  data = d,#rbind(d1, d2), 
  aes(x = x, y = y, fill = as.factor(attitude))
) +
  geom_rect(
    aes(xmin = 0.85, xmax = 1.45, ymin = -0.2, ymax = 1.7), 
    fill = "gray90"#"gray75"
  ) +
  geom_rect(
    aes(xmin = 1.55, xmax = 2.25, ymin = -0.2, ymax = 1.7), 
    fill = "gray90"
  ) +
  geom_text(
    check_overlap = TRUE,
    aes(x = 1.15, y = 1.5),
    label = "ethnic group A",
    size = 3
  ) +
  geom_text(
    check_overlap = TRUE,
    aes(x = 1.9, y = 1.5),
    label = "ethnic group B",
    size = 3
  ) +
  geom_jitter(shape = 21, size = 4, width = 0.04, height = 0.04) +
  facet_grid(~condition, scales = "free_x") +
  scale_y_continuous(limits = c(-0.5, 1.8), expand = c(0,0)) +
  scale_x_continuous(expand = c(0.1,0.1)) +
  scale_fill_manual(
    values = c("-" = "white", "+" = "gray20"), label=c("", "")) +
  labs(fill = "attitudes:") +
  theme(
    panel.background = element_rect(fill = "gray97"),
    plot.background = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    strip.background = element_rect(fill = "gray97"),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.key = element_blank(),
    axis.title= element_blank(),
    legend.title = element_text(size = 9),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
dev.off()

rm(d, d1, d2)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Chapter 2
#
# RI-model dominos 
worldSize = 10
iterations = 30
H = 0.6
mu = 0.5 # max opinion swing
seed = 222#111#12345#sample(-999999:999999, size = 1)


set.seed(seed)

# opinion matrix (dynamic)
m <- matrix(
  runif(min = 0, max = 1, n = worldSize ^ 2),
  nrow = worldSize,
  ncol = worldSize
)

# we can save the dynamics as an array:
w <- array(NA, dim = c(worldSize, worldSize, iterations + 1))
w[,,1] <- m

#group matrix (static)
gm <- matrix(1, nrow = worldSize, ncol = worldSize)
gm[,(round(worldSize / 2) + 1):worldSize] <- 0




print(paste("Running seed", seed))
pb <- txtProgressBar(min = 1, max = iterations, style = 3)
enableJIT(1)
for (t in 1:iterations) {
  setTxtProgressBar(pb, t)
  
  #shuffling order in which agents are called:
  shuffAgents <- arrayInd(
    ind = sample(x = 1:(worldSize ^ 2), size = worldSize ^ 2, replace = FALSE), 
    .dim = c(worldSize, worldSize)
  )
  
  # for all agents in random sequence:
  for (i in 1:(worldSize ^ 2)) {
    x <- shuffAgents[i,2]
    y <- shuffAgents[i,1]
    
    # find interaction partner j
    neighbors <- mooreNeigh(x = x, y = y, maxx = worldSize, maxy = worldSize)
    j <- neighbors[sample(1:nrow(neighbors), size = 1),]
    
    oi <- m[y,x]
    gi <- gm[y,x]
    
    oj <- m[j$y, j$x]
    gj <- gm[j$y, j$x]
    
    weight <- calcW(oi, oj, gi, gj, H = H, mu = mu)
    
    m[y,x] <- truncate(newO(oi = oi, oj = oj, w = weight), min = 0, max = 1 )
  }
  
  # Saving current time step
  w[,,t + 1] <- m
}
close(pb)
enableJIT(0)

image(t(w[,,1]))
image(t(w[,,iterations + 1]))


################################################################################

# take the last iteration stored in w:
d <- reshape2::melt(w[,,dim(w)[3]])
names(d) <- c("y", "x", "o")

#d <- subset(d, d$x %in% 1:size & d$y %in% 1:size)

ggplot(d, aes(x = x, y = y, fill = o)) +
  geom_tile() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_reverse(expand = c(0,0)) +
  #ggtitle("Average of rules 0-255") +
  scale_fill_viridis_c(option = "A") +
  theme(
    panel.background = element_blank(), 
    panel.grid = element_blank(),     
    axis.ticks = element_blank(), 
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(vjust = -1)
  )

################################################################################
# Plotting in 3d with rgl

for (i in 1:nrow(d)) ifelse(
  gm[d$y[i], d$x[i]],
  d$color[i] <- "gray80",
  d$color[i] <- "orange"#"#ffc569"#"orange"
)

printBox <- \(x, y, z, x1, y1, z1, o = 0.5, color = "gray40", strength = 2) {
  i <- cube3d(col = color, alpha = 1, shininess = 3) # create a cube as mesh object
  i <- scale3d(i, x1, y1, z1)   # now scale that object by x1,y1,z1
  i <- translate3d(i, x, y, z)  # now move it to x,y,z
  i <- rotate3d(
    i,
    angle = (pi/4 - (o * pi / 2)) * strength * -1,
    x = x,#x - 10,#-1,
    y = y,#y,#-2,
    z = -1 #0.3
  )
  shade3d(i)
  
}

height <- 1.5
side <- 0.3

open3d()

clear3d(type = "lights")
#light3d(x = 1, y = 1, z = 1, ambient = "red")
light3d(theta = -50, phi = 30, diffuse = "white", specular = "gray50")
#light3d(theta = -25, phi = 30, diffuse = "white", specular = "gray50")

#cube3d(col = color, alpha = 1, shininess = 3) # create a cube as mesh object
#i <- scale3d(i, x1, y1, z1)   # now scale that object by x1,y1,z1
#i <- translate3d(i, x, y, z)



for (i in 1:nrow(d)) {
  printBox(
    x = d$x[i], y = d$y[i], z = 0, # coordinates
    x1 = side, y1 = side, z1 = height, # dimension
    color = d$color[i], o = d$o[i], strength = 1 # aesthetics
  )
}

#bbox3d()
#rgl.bbox()
#axes3d()

snapshot3d("./output/Chapter_2.png", fmt = "png", width = 2000, height = 1500)
close3d()


rm(
  d, gm, j, m, neighbors, pb, shuffAgents, gi, gj, H, height, i, iterations,
  mu, oi, oj, seed, side, t, w, weight, worldSize, x, y, printBox
)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Chapter 3
#
# Network layout spelling out words


seed = 123456 #1234 #12345
set.seed(seed)

# Using ggplot to render text.
# The rendered string in a raster image will serve as our density maps to
# determine where the nodes will be.
png(
  filename = "./output/inputYes.png",
  width = 100, height = 60, res = 400, units = "px"
)
ggplot() +
  geom_text(aes(x = 0, y = 1, label = "yes")) +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(-4, -3, -3, -5), "pt"),
  )
dev.off()

png(
  filename = "./output/inputNo.png",
  width = 100, height = 60, res = 400, units = "px"
)
ggplot() +
  geom_text(aes(x = 0, y = 1, label = "no")) +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(-4, -3, -3, -5), "pt"),
  )
dev.off()

# Now we re-import these density maps as matrices.
yes <- png::readPNG("./output/inputYes.png")
yes <- 1 - ((yes[,,1] + yes[,,2] + yes[,,3]) / 3)
yes <- yes[60:1,]
no <- png::readPNG("./output/inputNo.png")
no <- 1 - ((no[,,1] + no[,,2] + no[,,3]) / 3)
no <- no[60:1,]

# Picking some random points in the plane. I use Latin hypercube sampling to
# ensure they're evenly spread out.
density = 700

lhs <- lhs::randomLHS(n = density, k = 2)
dYes <- data.frame(
  x = (lhs[,1] * 100) + 1,
  y = (lhs[,2] * 60) + 1,
  z = 0
)
lhs <- lhs::randomLHS(n = density, k = 2)
dNo <- data.frame(
  x = (lhs[,1] * 100) + 1,
  y = (lhs[,2] * 60) + 1,
  z = 0
)

# Removing the points that are not on the letters
for (i in 1:nrow(dYes)) {
  if(yes[floor(dYes$y[i]), floor(dYes$x[i])] != 0) dYes$z[i] <- 1
}
for (i in 1:nrow(dNo)) {
  if(no[floor(dNo$y[i]), floor(dNo$x[i])] != 0) dNo$z[i] <- 1
}
dYes <- dYes[dYes$z == 1, c("x", "y")]
dNo <- dNo[dNo$z == 1, c("x", "y")]
plot(dYes$x, dYes$y); plot(dNo$x, dNo$y)

# some of these points are too close to each other and clutter the plot.
# I get rid of them -- in a very inefficient way:
#dYes$remove <- FALSE
ddist <- as.matrix(dist(x = dYes))
diag(ddist) <- NA
howManyToRemove = round(nrow(dYes) * 0.55)
for (i in 1:howManyToRemove) {
  rem <- which(ddist == min(ddist, na.rm = TRUE), arr.ind = TRUE)[1,1]
  dYes <- dYes[-rem,]
  ddist <- as.matrix(dist(x = dYes))
  diag(ddist) <- NA
}

ddist <- as.matrix(dist(x = dNo))
diag(ddist) <- NA
howManyToRemove = round(nrow(dNo) * 0.55)
for (i in 1:howManyToRemove) {
  rem <- which(ddist == min(ddist, na.rm = TRUE), arr.ind = TRUE)[1,1]
  dNo <- dNo[-rem,]
  ddist <- as.matrix(dist(x = dNo))
  diag(ddist) <- NA
}


# Now we can make a network by creating ties between nearby nodes.
# Then we plot it.
proximityThreshold = 9

ddistY <- as.matrix(dist(x = dYes))
ddistN <- as.matrix(dist(x = dNo))

for(r in 1:nrow(ddistY)) for(c in 1:nrow(ddistY)) ifelse(
  ddistY[r,c] > proximityThreshold, ddistY[r,c] <- 0, ddistY[r,c] <- 1
)
for(r in 1:nrow(ddistN)) for(c in 1:nrow(ddistN)) ifelse(
  ddistN[r,c] > proximityThreshold, ddistN[r,c] <- 0, ddistN[r,c] <- 1
)

diag(ddistY) <- 0
diag(ddistN) <- 0

gY <- igraph::graph_from_adjacency_matrix(
  adjmatrix = ddistY,
  mode = "undirected"
)
gN <- igraph::graph_from_adjacency_matrix(
  adjmatrix = ddistN,
  mode = "undirected"
)
#plot(g)

png(
  filename = "./output/Chapter_3a.png",
  height = 750, width = 1200, res = 110, units = "px"
)
pY <- plot.igraph(
  gY,
  #rescale = FALSE,
  asp = 0.5,
  vertex.size = 5,
  vertex.color = "gray40", 
  vertex.frame.color = "gray30",
  edge.color = "gray50",
  edge.curved = FALSE, 
  vertex.label = NA,
  layout = as.matrix(dYes)
)
dev.off()

png(
  filename = "./output/Chapter_3b.png",
  height = 750, width = 900, res = 110, units = "px"
)
plot.igraph(
  gN,
  asp = 0.5,
  vertex.size = 5,
  vertex.color = "orange", 
  vertex.frame.color = "darkorange",
  edge.color = "gray50",
  edge.curved = FALSE, 
  vertex.label = NA,
  layout = as.matrix(dNo)
)
dev.off()

#ggpubr::ggarrange(pY, pN, ncol = 1)
rm(
  ddist, ddistN, ddistY, dNo, dYes, gN, gY, lhs, no, yes,
  c, density, howManyToRemove, i, proximityThreshold, pY, r, rem, seed
)

# We don't need the images anymore:
file.remove(paste0("./output/", c("inputYes", "inputNo"), ".png"))


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Chapter 4
#
# Segregated streets -- inspired by: https://github.com/marcusvolz/metropolis

set.seed(123)

createCity <- function( 
  n = 1000,
  r = 0.02,
  aNoise = 1,
  pBranch = 0.2,
  seed = sample(0:99999, size = 1)
) {
  set.seed(seed)
  
  points <- data.frame(
    x = numeric(n),
    y = numeric(n),
    dir = numeric(n),
    level = integer(n),
    root = integer(n)
  )
  edges <-  data.frame(
    x = numeric(n),
    y = numeric(n),
    xend = numeric(n),
    yend = numeric(n),
    level = integer(n),
    root = integer(n)
  )
  
  points$x[1:2] <- c(0.4, 0.6) #runif(2, 0, width)
  points$y[1:2] <- c(0.3, 0.3) #runif(2, 0, width)
  points$dir[1:2] <- runif(2, -2*pi, 2*pi)
  points$level[1:2] <- 1
  points$root[1:2] <- c(1, 2)
  
  cat("City generation...\n")
  progressBar <- txtProgressBar(min = 0, max = n, style = 3)
  enableJIT(1)
  for (i in (3):n) {
    attempts <- 0
    while (attempts <= 50) {
      attempts <- attempts + 1
      random_point <- points[sample(1:(i - 1), size = 1),]
      branch <- sample(
        c(-1, 0, 1), size = 1,
        prob = c(pBranch / 2, 1 - pBranch, pBranch / 2)
      )
      
      delta <- aNoise * 2 * pi / 180
      alpha <- random_point$dir + runif(1, -(delta), delta) + (branch * pi/2)
      ifelse(
        branch != 0,
        newLevel <- random_point$level + 1,
        newLevel <- random_point$level
      )
      v <- c(cos(alpha), sin(alpha)) * r * (1 + (1 / newLevel))
      xj <- random_point$x + v[1]
      yj <- random_point$y + v[2]
      
      if(xj < 0 | xj > 1 | yj < 0 | yj > 1) next
      
      distances <- calcDist(
        xj, points$x[1:i], yj, points$y[1:i], fastMode = FALSE)
      if (min(distances) >= r) {
        points[i, ] <- c(xj, yj, alpha, newLevel, random_point$root)
        edges[i, ] <- 
          c(xj, yj, random_point$x, random_point$y, newLevel, random_point$root)
      }
    }
    setTxtProgressBar(progressBar, i)
  }
  enableJIT(0)
  
  edges <- edges[edges$level != 0 & edges$root != 0,]
  points <- points[points$root != 0,]
  
  edges$faultline <- edges$meany <- edges$meanx <- NA
  edges$meanx <- (edges$x + edges$xend) / 2
  edges$meany <- (edges$y + edges$yend) / 2

  
  a <- which(edges$root == 1)
  b <- which(edges$root != 1)
  for (i in a) {
    edges$faultline[i] <- 1 - min(calcDist(
      edges$meanx[i], edges$meanx[b],
      edges$meany[i], edges$meany[b],
      fastMode = FALSE
    ))
  }
  for (i in b) {
    edges$faultline[i] <- min(calcDist(
      edges$meanx[i], edges$meanx[a],
      edges$meany[i], edges$meany[a],
      fastMode = FALSE
    ))
  }
  #plot(edges$meanx, edges$faultline)
  
  close(progressBar)
  return(list(
    edges = edges,
    points = points,
    parameters = list(
      n = as.integer(n),
      r = r,
      pBranch = pBranch,
      seed = seed
    )
  ))
}


plotCity <- function(edges = edges, streetWeight = 1.1) {
  p <- ggplot(edges, aes(x = x, y = y)) +
    geom_segment( # Background
      aes(xend = xend, yend = yend),
      color = "gray50", size = streetWeight, lineend = "round"
    ) +
    geom_segment( # Showing faultlines
      aes(xend = xend, yend = yend, color = faultline),
      size = streetWeight * 0.9, lineend = "round"
    ) +
    scale_color_gradientn(
      colors = c(
        "orange",# "darkorange",
        "gray90", "gray90", "gray90",
        "black"#, "black"
      )
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    coord_equal() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
  return(p)
}

#_______________________________________________________________________________


c <- createCity(n = 5000, r = 0.02, aNoise = 2, pBranch = 0.05)
list2env(c, globalenv())
rm(c)


print(plotCity(edges, streetWeight = 1.1))
print(parameters$seed)


png(
  filename = "./output/Chapter_4.png",
  width = 2000, height = 2000, res = 150, units = "px"
)
print(plotCity(edges, streetWeight = 4))
dev.off()

rm(edges, parameters, points, plotCity, createCity)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Chapter 5
#
# Parameter space

set.seed(123)


d <- ambient::noise_perlin(
  dim = c(10, 10, 10),
  frequency = 0.05,
  gain = 0.5,
  pertubation = "fractal"
) |> reshape2::melt()

names(d) <- c("x", "y", "z", "color")

# normalize and discretize the color values:
d$color <- (d$color - min(d$color)) / (max(d$color) - min(d$color))
#d$color <- findInterval(x = d$color, vec = 1:3/4) 
for (i in 1:nrow(d)) ifelse(d$color[i] > 0.6, d$color[i] <- 0, d$color[i] <- 1)

ax <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)


p <- plotly::plot_ly(d) |>
  plotly::add_markers(
    x = ~x,
    y = ~y,
    z = ~z, 
    color = ~color,
    type = "scatter3d",
    size = I(3300),#I(2100),
    marker = list(opacity = 0.5, line = list(width = 0)),
    #colors = viridis::viridis_pal(option = "B", begin = 0.6, end = 0.85) (3)
    #colors = viridis::viridis_pal(option = "B", begin = 0.65, end = 0.9) (3)
    colors = c("orange", "red", "gray70")
  ) |>
  plotly::layout(scene = list(
    xaxis = ax,yaxis = ax, zaxis = ax,
    #camera = list(eye = list(x = 1.8, y = 1.2, z = -0.3))
    camera = list(eye = list(x = 1.5, y = 2, z = 1.7))
  )) |>
  plotly::hide_colorbar()
print(p)

#plotly::kaleido(
plotly::orca(
  p, 
  file = "./output/Chapter_5.png",
  width = 2000,
  height = 2000
)
#plotly::save_image(p, file = "./output/parameter_space.png")

rm(ax, d, p, i)
