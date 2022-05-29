install.packages("networkD3")
library(networkD3)
install.packages("ggalluvial")
library(ggalluvial)

Pat_data3 = read.delim("patdata_fungi_sneak.txt",sep='\t', header = T, row.names = 1)
my_data3 = Pat_data3


is_alluvia_form(as.data.frame(Pat_data3), axes = 1:3, silent = TRUE)


ggplot(as.data.frame(my_data3),
       aes(axis1 = age_group, axis2 = Culture_fungi)) +
  geom_alluvium(aes(fill = ICU_discharge_alive), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("age_group", "Culture_fungi"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("ICU discharge status by sex, age, and cultural status")


ggplot(as.data.frame(my_data3),
       aes(axis1 = sex, axis2 = age_group, axis3 = Culture_fungi)) +
  geom_alluvium(aes(fill = ICU_discharge_alive), width = 1/12) +
  geom_stratum(width = 1/12, fill = "orange", color = "gray") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sex", "age_group", "Culture_fungi"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("ICU discharge status by sex, age, and cultural status")


ggplot(as.data.frame(my_data3),
       aes(axis1 = sex, axis2 = Leucocytes_range, axis3 = CRP_range, axis4 = PCT_range)) +
  geom_alluvium(aes(fill = ICU_discharge_alive), width = 1/12) +
  geom_stratum(width = 1/12, fill = "orange", color = "gray") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sex", "Leucocytes_range", "CRP_range", "PCT_range"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("ICU discharge status by sex, Leucocytes, CRP and PCT levels")


ggplot(as.data.frame(my_data3),
       aes(axis1 = Leucocytes_range, axis2 = CRP_range, axis3 = PCT_range)) +
  geom_alluvium(aes(fill = ICU_discharge_alive), width = 1/12) +
  geom_stratum(width = 1/12, fill = "orange", color = "gray") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Leucocytes_range", "CRP_range", "PCT_range"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("ICU discharge status by Leucocytes, CRP and PCT levels")


ggplot(as.data.frame(my_data3),
       aes(axis1 = Leucocytes_range, axis2 = Micro_Granulocytes, axis3 = Culture_fungi)) +
  geom_alluvium(aes(fill = ICU_discharge_alive), width = 1/12) +
  geom_stratum(width = 1/12, fill = "orange", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Leucocytes_range", "Micro_Granulocytes", "Culture_fungi"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("ICU discharge status by Leucocytes, Microscopic Granulocytes, and cultural status")





ggplot(as.data.frame(my_data3),
       aes(axis1 = sex, axis2 = age_group, axis3 = Six.day_evaluation)) +
  geom_alluvium(aes(fill = ICU_discharge_alive), width = 1/12) +
  geom_stratum(width = 1/12, fill = "orange", color = "black") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sex", "age_group", "6day evaluation"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("ICU discharge status by sex, age, and Six day evaluation")
 
 ####################################################
various statistics
#####################################################

ggbetweenstats(
data  = my_data3,
x     = Culture_fungi,
y     = Six.day_evaluation,
type = "n",
title = "Six day evaluation of fungal culture positive and negative samples"
)

grouped_ggbetweenstats(
  data             = my_data3,
  x     = age_group,
  y     = Six.day_evaluation,
  grouping.var     = sex, ## grouping variable
  outlier.tagging  = TRUE, ## whether outliers need to be tagged
  outlier.coef     = 2,
  ggsignif.args    = list(textsize = 4, tip_length = 0.01),
  p.adjust.method  = "bonferroni", ## method for adjusting p-values for multiple comparisons
  ## adding new components to `{ggstatsplot}` default
  ggplot.component = list(ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis())),
  caption          = "INTeGRATE",
  palette          = "default_jama",
  package          = "ggsci",
  plotgrid.args    = list(nrow = 1),
  annotation.args  = list(title = "Differences in movie length by mpaa ratings for different genres")
)

ggbetweenstats(
  data  = my_data3,
  x     = age_group,
  y     = Six.day_evaluation,
  type = "n",
  title = "Six day evaluation of different age groups"
)

ggbetweenstats(
  data  = my_data3,
  x     = Antifungal_thrapy.5d.,
  y     = Six.day_evaluation,
  type = "n",
  title = "Six day evaluation in antifungal treatment groups"
)


ggpiestats(
  data         = mtcars,
  x            = am,
  y            = cyl,
  package      = "wesanderson",
  palette      = "Royal1",
  title        = "Dataset: Motor Trend Car Road Tests", ## title for the plot
  legend.title = "Transmission", ## title for the legend
  caption      = "Source: 1974 Motor Trend US magazine"
)




ggbarstats(
  data             = my_data3,
  x                = sex,
  y                = ICU_discharge_alive,
  title            = "ICU discharge by sex",
  xlab             = "ICU discharge alive",
  legend.title     = "sex",
  ggtheme          = hrbrthemes::theme_ipsum_pub(),
  ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
  palette          = "Set2"
)



grouped_ggbarstats(
  data         = my_data3,
  x            = ICU_discharge_alive,
  y            = age_group,
  grouping.var = sex,
  package      = "wesanderson",
  palette      = "Darjeeling2",
  ggtheme      = ggthemes::theme_tufte(base_size = 12)
)
