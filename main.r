mu_0 <- 4 * pi * 10^-7  
epsilon_0 <- 8.854 * 10^-12  
c <- 299792458

calc_electric_field <- function(charge, position, observation_point) {
  r <- observation_point - position
  norm_r <- sqrt(sum(r^2))
  
  k <- 1 / (4 * pi * epsilon_0)
  E <- (k * charge * r) / (norm_r^3)
  
  E
}

calc_magnetic_field <- function(current, position, observation_point) {
  r <- observation_point - position
  norm_r <- sqrt(sum(r^2))
  
  k <- mu_0 / (4 * pi)
  cross <- crossprod(current, r)
  B <- k * cross / (norm_r^2)
  
  B
}

calc_force <- function(charge1, position1, charge2, position2) {
  r <- position2 - position1
  norm_r <- sqrt(sum(r^2))
  
  k <- 1 / (4 * pi * epsilon_0)
  F <- (k * charge1 * charge2 * r) / (norm_r^3)
  
  F
}

calc_electric_field_distribution <- function(charge, position, observation_point) {
  n <- nrow(charge)
  E_total <- c(0, 0, 0)
  
  for (i in 1:n) {
    E_total <- E_total + calc_electric_field(charge[i], position[i,], observation_point)
  }
  
  E_total
}

calc_magnetic_field_distribution <- function(current, position, observation_point) {
  n <- nrow(current)
  B_total <- c(0, 0, 0)
  
  for (i in 1:n) {
    B_total <- B_total + calc_magnetic_field(current[i,], position[i,], observation_point)
  }
  
  B_total
}

calc_wire_force <- function(current1, position1, current2, position2, length) {
  B1 <- calc_magnetic_field(current1, position1, position2)
  B2 <- calc_magnetic_field(current2, position2, position1)
  
  F <- length * crossprod(B1, B2) / (4 * pi * mu_0)
  
  F
}

calc_resistance <- function(length, area, conductivity) {
  R <- length / (area * conductivity)
  
  R
}

calc_capacitance <- function(area, distance) {
  C <- epsilon_0 * area / distance
  
  C
}

calc_inductance <- function(n_turns, length, radius) {
  L <- mu_0 * n_turns^2 * pi * radius^2 * length
  
  L
}
library(rgl)

calculate_electric_field <- function(r, t) {
  c <- 299792458 
  q <- 1.602176634e-19  
  eps0 <- 8.8541878128e-12 
  omega <- 2 * pi * c / 800e-9
  
  dipole_position <- c(0, 0, 0)
  dipole_moment <- c(0, 0, 1e-29)
  distance <- sqrt(sum((r - dipole_position)^2))
  
  unit_vector <- (r - dipole_position) / distance
  
  amplitude <- (q / (4 * pi * eps0 * distance^3)) * sqrt(3 / 4) * omega * abs(dipole_moment) * sin(omega * t)
  electric_field <- amplitude * unit_vector
  
  electric_field
}

calculate_intensity <- function(r, t) {
  electric_field <- calculate_electric_field(r, t)
  intensity <- sum(electric_field^2) / (2 * 376.73031346177)
  intensity
}

calculate_intensity_map <- function(x_range, y_range, z_range, t) {
  num_x <- length(x_range)
  num_y <- length(y_range)
  num_z <- length(z_range)
  intensity_map <- array(0, dim = c(num_x, num_y, num_z))
  
  for (i in 1:num_x) {
    for (j in 1:num_y) {
      for (k in 1:num_z) {
        r <- c(x_range[i], y_range[j], z_range[k])
        intensity_map[i, j, k] <- calculate_intensity(r, t)
      }
    }
  }
  
  intensity_map
}

plot_intensity_map <- function(intensity_map) {
  num_x <- dim(intensity_map)[1]
  num_y <- dim(intensity_map)[2]
  num_z <- dim(intensity_map)[3]
  max_intensity <- max(intensity_map)
  
  for (i in 1:num_x) {
    for (j in 1:num_y) {
      for (k in 1:num_z) {
        r <- c(i, j, k)
        intensity <- intensity_map[i, j, k]
        
        col <- heat.colors(256)[as.integer(intensity / max_intensity * 255) + 1]
        points3D(r, col = col, add = TRUE, pch = 16)
      }
    }
  }
  
  title(main = "Intensity Map", xlab = "X", ylab = "Y", zlab = "Z")
}

animate_intensity_map <- function(x_range, y_range, z_range, num_frames) {
  max_time <- 5e-15
  times <- seq(0, max_time, length.out = num_frames)
  intensity_maps <- list()
  
  for (i in 1:length(times)) {
    t <- times[i]
    intensity_map <- calculate_intensity_map(x_range, y_range, z_range, t)
    intensity_maps[[i]] <- intensity_map
  }
  
  animation::ani.options(interval = 1/30)
  saveGIF({
    for (i in 1:length(intensity_maps)) {
      intensity_map <- intensity_maps[[i]]
      plot_intensity_map(intensity_map)
    }
  }, movie.name = "intensity_map.gif", ani.width = 600, ani.height = 600)
}