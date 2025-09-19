library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

list.files(path = ".", full.names = TRUE) %>%
  purrr::walk(~file.rename(.x, sub("Placeholder 1", "", .x)))
