#!/usr/bin/env Rscript
# Script to generate pipeML hex sticker logo

library(hexSticker)
library(ggplot2)
library(showtext)
library(sysfonts)
library(grid)

# Add a nice Google font
font_add_google("Fira Sans", "firasans")
font_add_google("Roboto Mono", "robotomono")
showtext_auto()

# ---- Build the inner plot: gear + circuit/ML nodes ----

# Gear shape (outer gear with teeth)
gear <- function(cx, cy, r_inner, r_outer, n_teeth = 8) {
  angles <- seq(0, 2 * pi, length.out = n_teeth * 4 + 1)
  r <- rep(c(r_outer, r_outer, r_inner, r_inner), n_teeth)
  r <- c(r, r[1])
  data.frame(
    x = cx + r * cos(angles),
    y = cy + r * sin(angles)
  )
}

# Create gear polygon
gear_df <- gear(0, 0, r_inner = 0.65, r_outer = 0.85, n_teeth = 8)

# Inner circle for the gear hub
circle_points <- function(cx, cy, r, n = 100) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

inner_circle <- circle_points(0, 0, 0.4)
hub_circle <- circle_points(0, 0, 0.15)

# Neural network / circuit nodes
nodes <- data.frame(
  x = c(-0.25, 0.25, 0, -0.2, 0.2),
  y = c(0.15, 0.15, -0.15, -0.15, -0.15)
)

# Connections between nodes (simple network)
edges <- data.frame(
  x    = c(-0.25, -0.25, 0.25, 0.25, 0, 0),
  y    = c(0.15,  0.15,  0.15, 0.15, -0.15, -0.15),
  xend = c(0,    -0.2,   0,    0.2,  -0.2,  0.2),
  yend = c(-0.15, -0.15, -0.15, -0.15, -0.15, -0.15)
)

# Pipeline arrows running left-to-right through the gear
pipe_df <- data.frame(
  x    = c(-1.1, -0.15),
  y    = c(0, 0),
  xend = c(-0.15, 1.1),
  yend = c(0, 0)
)

# Small circuit branch dots
circuit_dots <- data.frame(
  x = c(0.9, 0.75, 1.0),
  y = c(0.25, -0.25, -0.15)
)

circuit_lines <- data.frame(
  x    = c(0.7, 0.7, 0.85),
  y    = c(0, 0, 0),
  xend = c(0.9, 0.75, 1.0),
  yend = c(0.25, -0.25, -0.15)
)

# Build the subplot
p <- ggplot() +
  # Gear body
  geom_polygon(data = gear_df, aes(x, y),
               fill = "#3A7CA5", color = "#E8E8E8", linewidth = 0.8) +
  # Inner circle (gear cutout)
  geom_polygon(data = inner_circle, aes(x, y),
               fill = "#1B2A4A", color = "#E8E8E8", linewidth = 0.5) +
  # Pipeline arrows
  geom_segment(data = pipe_df, aes(x = x, y = y, xend = xend, yend = yend),
               color = "#E8E8E8", linewidth = 1.2,
               arrow = arrow(length = unit(0.06, "npc"), type = "closed")) +
  # Neural network edges inside gear
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
               color = "#81D4FA", linewidth = 0.6) +
  # Neural network nodes inside gear
  geom_point(data = nodes, aes(x, y),
             color = "#E8E8E8", fill = "#FF7043", shape = 21, size = 3, stroke = 0.5) +
  # Hub
  geom_polygon(data = hub_circle, aes(x, y),
               fill = "#FF7043", color = "#E8E8E8", linewidth = 0.4) +
  # Circuit branches out of pipe
  geom_segment(data = circuit_lines, aes(x = x, y = y, xend = xend, yend = yend),
               color = "#81D4FA", linewidth = 0.5) +
  geom_point(data = circuit_dots, aes(x, y),
             color = "#81D4FA", size = 1.5) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.0, 1.0)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# ---- Create hex sticker ----
sticker(
  subplot = p,
  package = "pipeML",
  # Package name styling
  p_size = 24,
  p_color = "#E8E8E8",
  p_family = "firasans",
  p_fontface = "bold",
  p_y = 1.45,
  # Subplot placement
  s_x = 1, s_y = 0.78,
  s_width = 1.6, s_height = 1.2,
  # Hex styling
  h_fill = "#1B2A4A",
  h_color = "#3A7CA5",
  h_size = 1.8,
  # White border
  white_around_sticker = FALSE,
  # Output
  filename = "man/figures/logo.png",
  dpi = 300
)

cat("Hex sticker created at man/figures/logo.png\n")
