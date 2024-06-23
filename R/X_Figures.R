
# X_Figures ####

library(tidyverse); library(cowplot); library(fs); library(magrittr)

theme_set(theme_cowplot())

dir_create("Figures")

# Figure 1: Epidemiological Simulations ####

# Panel A: Growth Curves ####

TimeStepDF <- readRDS("Data/Outputs/PopulationTimeSteps.rds")

(TimeFigure <- 
   TimeStepDF %>%
   # mutate(Population = paste0("Group: ", Population)) %>% 
   # mutate(Population = paste0("Group: ", Population)) %>% 
   mutate(File = paste(Rep, R, P_I, sep = "_")) %>% 
   mutate_at("Year", ~factor(.x, levels = c(2013:2017, "\U2022 Maria \U2022", 2018, 2019, 2021, 2022))) %>% 
   ggplot(aes(Time + 1, PropInf, colour = Year)) + 
   geom_line(aes(group = File), alpha = 0.05) + 
   geom_line(data = . %>% filter(Time > 10000), aes(group = File), alpha = 1) + 
   scale_y_continuous(breaks = c(0:2/2)) +
   scale_x_continuous(limits = c(0, 1500)) +
   scale_colour_discrete_divergingx(palette = "PuOR", 
                                    limits = c(2013:2017, "\U2022 Maria \U2022", 2018, 2019, 2021, 2022)) +
   labs(colour = NULL, x = "Time step", y = "Proportion infected") +
   theme(legend.position = "top", legend.justification = 0.5) +
   facet_grid(Population ~ .) + 
   theme(strip.background = element_blank()) +#,
   guides(colour = guide_legend(reverse = F,
                                direction = "horizontal",
                                label.position = "top",
                                label.theme = element_text(angle = 0), nrow = 1)) +
   NULL)

# Panel B: Camille Violins ####


# Panel C: P_I versus mean infection timestep ####

OutputDF <- readRDS("Data/Outputs/PopulationTimes.rds")

OutputDF %<>% 
  mutate_at("Rep", str_trim) %>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate_at("Rep", ~str_replace(.x, "HH", "H")) %>% 
  mutate_at("Rep", ~str_replace(.x, "TT", "T")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5))

OutputDF %<>% mutate(PostMaria = as.factor(as.numeric(Year %in% c(2018:2021))))

OutputDF %<>% mutate_at(c("Mean", "Max"), ~log(.x + 1))

OutputDF %<>% filter(between(Mean, 3, 7))

(P_I_Figure <- 
    OutputDF %>% 
    mutate_at("Year", ~factor(.x, levels = c(2013:2017, "\U2022 Maria \U2022", 2018, 2019, 2021, 2022))) %>% 
    ggplot(aes(P_I, exp(Mean), colour = Year, group = Rep)) + 
    geom_point(alpha = 0.2) +
    geom_smooth(method = lm, fill = NA) +
    scale_y_log10(limits = c(10, NA)) +
    labs(y = "Mean infection timestep", x = "Pathogen infectivity") +
    scale_colour_discrete_divergingx(palette = "PuOR", 
                                     limits = c(2013:2017, "[ Maria ]", 2018, 2019, 2021, 2022)) +    
    theme(legend.position = "none"))

# Panel D: Model effect comparisons ####

MCMC1 <- readRDS("Output/EpidemiologyModel.rds")
MCMC1 <- readRDS("Model.rds")

MCMCOutput <- 
  MCMC1$Sol %>% data.frame %>% 
  mutate_at("P_I", ~.x/10) %>% 
  dplyr::select(P_I, PostMaria1) %>% 
  pivot_longer(1:2, names_to = "rowname", values_to = "Estimate")

CoefDF <- MCMC1 %>% summary %>% extract2(5) %>% as.data.frame %>% 
  rename(Lower = 2, Upper = 3, Estimate = 1) %>% 
  slice(2:3)

CoefDF[1, c("Estimate", "Lower", "Upper")] <-
  CoefDF[1, c("Estimate", "Lower", "Upper")]/10

(EffectPlot <- 
    CoefDF %>% 
    rownames_to_column() %>% 
    ggplot(aes(rowname, Estimate)) +
    geom_violin(data = MCMCOutput, scale = "width", fill = AlberColours[[1]], alpha = 0.3, colour = NA) +
    geom_hline(lty = 2, alpha = 0.3, yintercept = 0) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5) +
    geom_point(size = 2) +
    geom_point(colour = "white") +
    scale_x_discrete(labels = c("+10% Infectivity", "Hurricane")) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(y = "Effect estimate", x = NULL))

# Combining ####

(TimeFigure)/
  ((P_I_Figure + 
      theme(axis.title.x = element_text(vjust = 10)) + 
      EffectPlot) + 
     plot_layout(widths = c(4, 1))) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure1.jpeg", units = "mm", 
       width = 250, height = 250, dpi = 600)

(TimeFigure)/
  ((P_I_Figure + 
      theme(axis.title.x = element_text(vjust = 10)) + 
      EffectPlot) + 
     plot_layout(widths = c(2, 1))) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure1Thin.jpeg", units = "mm", 
       width = 180, height = 250, dpi = 600)

ggsave("Figures/Figure1Thin.pdf", units = "mm", 
       width = 180, height = 250, dpi = 600)



# Old Figure Stuff #####

EdgeListList <- readRDS("Greg Data/TimeEdges.rds")

EdgeListList %<>%
  separate(DateTime, sep = " ", into = c("Date", "Time"))

EdgeListList %<>% mutate_at("Date", ~lubridate::ymd(.x))

EdgeListList %<>% 
  mutate(Year = str_split(Date, "-") %>% map_chr(1))

EdgeListList %>% 
  count(Rep, Date, Year) %>% 
  ggplot(aes(Rep, Date)) + 
  geom_tile(aes(fill = n)) + 
  coord_flip() +
  facet_wrap(~Year, scales = "free")

DiffDF <- 
  EdgeListList %>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  group_by(Rep, Date, Year) %>% 
  count %>% group_by(Rep, Year) %>% 
  mutate(DateGap = c(0, diff(as.numeric(Date))))

DiffDF %>% 
  filter(DateGap > 0) %>% 
  # mutate_at("Gap", ~(.x + 1)) %>% 
  SinaGraph("Rep", "DateGap", Just = T) + 
  scale_y_log10()

LongMaxes %>% 
  group_by(Rep) %>% 
  summarise_at(c("Mean", "Max"), mean) %>% 
  left_join(DiffDF %>% 
              group_by(Rep) %>% 
              summarise_at("DateGap", mean)) %>% 
  filter(DateGap < 4) %>% 
  ggplot(aes(DateGap, Mean)) +
  # geom_point() + 
  geom_text(aes(label = Rep)) +
  geom_smooth(method = lm)

# Plotting temporal windows ####

Window <- 30

SubEdgeList <-
  EdgeListList

SubEdgeList %<>% 
  group_by(Rep) %>% 
  mutate(NDate = Date - min(Date) + 1)

SubEdgeList %>% 
  ggplot(aes(Rep, NDate)) + 
  geom_point() + 
  coord_flip()

FullDays <- 
  max(SubEdgeList$NDate) - 
  min(SubEdgeList$NDate) - 
  Window

library(tidygraph)

BlankIndivs <- 
  SubEdgeList[,c("From", "To")] %>% unlist %>% unique %>% 
  data.frame(ID = .)

GraphList <- 
  0:FullDays %>% 
  map(function(Day){
    
    SubSubEdgeList <- 
      SubEdgeList %>% 
      filter(Rep == FocalRep) %>% 
      filter(NDate %in% (1:Window + Day)) %>% 
      dplyr::select(From, To) %>% 
      # graph_from_data_frame %>% #as_tbl_graph()
      tbl_graph(edges = ., 
                nodes = BlankIndivs) %>% 
      get.adjacency(sparse = F)
    
  })

library(ggregplot); library(ggraph); library(patchwork)

a <- GraphList[[1]]

Samples <- 
  sample(1:length(GraphList), 10) %>% 
  sort

GraphList[Samples] %>% 
  map(function(a){
    
    a %>% 
      ggraph("stress") + 
      geom_edge_link0() +
      geom_node_point() +
      coord_fixed()
    
  }) %>% ArrangeCowplot()
