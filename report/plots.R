library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

strong <- read.csv("~/OneDrive/github/KD-tree/time/strong/strong_v2.csv", header=TRUE)
weak <- read.csv("~/OneDrive/github/KD-tree/time/weak/weak.csv", header=TRUE)


plot_data <- ggplot(data = strong)+
  xlab("N process") +
  ylab("Time (s)") +
  theme_bw() + 
  # scale_colour_brewer(type = "seq", palette = "Set1") +
  scale_color_manual(values = c("#003f5c", "#bc5090", "#ffa600")) +
  scale_y_continuous(breaks = seq(0,650,by=20))+
  scale_x_continuous(trans="log2")+
  geom_line(mapping = aes(nprocs, serial, col = "Serial"), size = 1.5) +
  geom_point(mapping = aes(nprocs, serial, col = "Serial"), size = 3) +
  geom_line(mapping = aes(nprocs, mpi, col = "Open MPI"), size = 1.5) +
  geom_point(mapping = aes(nprocs, mpi, col = "Open MPI"), size = 3) +
  geom_line(mapping = aes(nprocs, omp, col = "OpenMP"), size = 1.5) +
  geom_point(mapping = aes(nprocs, omp, col = "OpenMP"), size = 3) + 
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = "top")
ggsave(paste0("~/OneDrive/github/KD-tree/report/strong.png"), width = 12, height = 6, dpi = 200)

# OpenMP

gg1 <- ggplot(data = strong)+
  xlab("N of process") +
  ylab("Time (s)") +
  theme_bw() + 
  scale_color_manual(values = "#bc5090") +
  geom_line(mapping = aes(nprocs, omp, col = "OpenMP"), size = 1.5) +
  geom_point(mapping = aes(nprocs, omp, col = "OpenMP"), size = 3) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = "top")
ggsave(paste0("~/OneDrive/github/KD-tree/report/strong_OMP.png"), width = 8, height = 6, dpi = 200)

gg2 <- ggplot(data = strong)+
  xlab("N of process") +
  ylab("Time (s)") +
  theme_bw() + 
  scale_color_manual(values = "#003f5c") +
  geom_line(mapping = aes(nprocs, mpi, col = "OpenMPI"), size = 1.5) +
  geom_point(mapping = aes(nprocs, mpi, col = "OpenMPI"), size = 3) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = "top")
ggsave(paste0("~/OneDrive/github/KD-tree/report/strong_MPI.png"), width = 8, height = 6, dpi = 200)

gg3 <- ggplot(data = strong)+
  xlab("N of process") +
  ylab("T(1)/T(N)") +
  theme_bw() + 
  scale_color_manual(values = c("#003f5c", "#bc5090")) +
  geom_line(mapping = aes(nprocs, speedupOMP, col = "Open MP"), size = 1.5) +
  geom_point(mapping = aes(nprocs, speedupOMP, col = "Open MP"), size = 3) +
  geom_line(mapping = aes(nprocs, speedupMPI, col = "OpenMPI"), size = 1.5) +
  geom_point(mapping = aes(nprocs, speedupMPI, col = "OpenMPI"), size = 3) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = "top")
ggsave(paste0("~/OneDrive/github/KD-tree/report/speedupStrong.png"), width = 12, height = 6, dpi = 200)

grid.arrange(gg1, gg2, ncol=2)

ggw1 <- ggplot(data = weak) +
  xlab("N of process (log)") +
  ylab("Time (s)") +
  theme_bw() + 
  scale_color_manual(values = c("#003f5c", "#bc5090")) +
  geom_line(mapping = aes(nprocs, mpi, col = "Open MPI"), size = 1.5) +
  geom_point(mapping = aes(nprocs, mpi, col = "Open MPI"), size = 3) +
  geom_line(mapping = aes(nprocs, omp, col = "OpenMP"), size = 1.5) +
  geom_point(mapping = aes(nprocs, omp, col = "OpenMP"), size = 3) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = "top")
ggsave(paste0("~/OneDrive/github/KD-tree/report/weak.png"), width = 12, height = 6, dpi = 200)

ggw2 <- ggplot(data = weak)+
  xlab("N of process (log)") +
  ylab("T(1)/T(N)") +
  theme_bw() + 
  scale_x_continuous(trans='log2') +
  scale_color_manual(values = c("#003f5c", "#bc5090")) +
  #scale_x_continuous(trans="log2")+
  geom_line(mapping = aes(nprocs, speedupMPI, col = "Open MPI"), size = 1.5) +
  geom_point(mapping = aes(nprocs, speedupMPI, col = "Open MPI"), size = 3) +
  geom_line(mapping = aes(nprocs, speedupOMP, col = "OpenMP"), size = 1.5) +
  geom_point(mapping = aes(nprocs, speedupOMP, col = "OpenMP"), size = 3) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(text = element_text(size = 20)) +
  theme(legend.position = "top") + 
  ylim(0, 5)
ggsave(paste0("~/OneDrive/github/KD-tree/report/speedupWeak.png"), width = 12, height = 6, dpi = 200)

grid.arrange(ggw1, ggw2, ncol=2)

