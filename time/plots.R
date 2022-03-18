library(dplyr)
library(ggplot2)
library(tidyr)

strong_mpi <- tibble(read.csv("~/OneDrive/github/KD-tree/time/strong/mpi.csv"))
strong_omp <- tibble(read.csv("~/OneDrive/github/KD-tree/time/strong/omp.csv"))
strong_serial <- tibble(read.csv("~/OneDrive/github/KD-tree/time/strong/serial.csv"))

colnames(strong_mpi) <- c("nprocs", "time_mpi")
colnames(strong_omp) <- c("nprocs", "time_omp")
colnames(strong_serial) <- c("nprocs", "time_serial")

dataset_scalability <- tibble(strong_mpi$nprocs, strong_mpi$time_mpi, 
                              strong_omp$time_omp, strong_serial$time_serial)
colnames(dataset_scalability) <- c("nprocs", "mpi", "omp", "serial")


plot_data <- ggplot(data = dataset_scalability)+
  xlab("N process") +
  ylab("Time (s)") +
  theme_bw() + 
  scale_colour_brewer(type = "seq", palette = "Set1") +
  scale_y_continuous(breaks = seq(0,650,by=20))+
  scale_x_continuous(trans="log2")+
  geom_line(mapping = aes(nprocs, serial, col = "Serial")) +
  geom_point(mapping = aes(nprocs, serial, col = "Serial")) +
  geom_line(mapping = aes(nprocs, mpi, col = "Open MPI")) +
  geom_point(mapping = aes(nprocs, mpi, col = "Open MPI")) +
  geom_line(mapping = aes(nprocs, omp, col = "OpenMP")) +
  geom_point(mapping = aes(nprocs, omp, col = "OpenMP"))
ggsave(paste0("~/OneDrive/github/KD-tree/time/strong.png"), width = 8, height = 6, dpi = 200)


