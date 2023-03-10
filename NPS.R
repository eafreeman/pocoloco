library(readxl)
library(tidyverse)

np <- read_excel("/Users/elizabethfreeman/Desktop/np_inj.xlsx")

np_long <- pivot_longer(np, cols = c(`day 0`, `day 2`, `day 4`, `day 7`), 
                        names_to = "day", values_to = "survived")

ggplot(data = np_long) +
  geom_line(aes(x = day, y = survived, group = concentration, color = concentration)) +
  facet_wrap(~sex)
  