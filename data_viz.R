
# Prepare data by linking ids to join genes/proteins expression readouts

library(dplyr)

data<- left_join(A, B, by="blockid")
data<- left_join(data, C, by="blockid")
data<- left_join(data, D, by="blockid")
rownames(data) <- data$blockid
data$blockid <- NULL

# Visualizations

# Upset plot for binary readouts

library(ComplexUpset)
library(ggplot2)
library(dplyr)

# Prepare upset_data

upset_data <- data

# Change readouts to binary values depending on agreed cutoffs
upset_data$A <- ifelse(upset_data$A >=25,25,upset_data$A)
upset_data <- upset_data %>% mutate (A=ifelse(upset_data$A==25,1,0))


upset_data$C <- ifelse(upset_data$C >=10,10,upset_data$C)
upset_data <- upset_data %>% mutate (C=ifelse(upset_data$C==10,1,0))

upset_data$D <- ifelse(upset_data$D >=10,10,upset_data$D)
upset_data <- upset_data %>% mutate (D=ifelse(upset_data$D==10,1,0))

upset_data$B <- ifelse(upset_data$B >=50,50,upset_data$B)
upset_data <- upset_data %>% mutate(B=ifelse(upset_data$B==50,1,0))

upset_data <- na.omit(upset_data)

size = get_size_mode('exclusive_intersection')
upset(
  upset_data,
  c("A", "B", "C", "D" ),
  name= "Marker Positivity",
  queries=list(
    upset_query(set='A', fill='#56B4E9'),
    upset_query(set='B', fill='#Ffb347') ,
    upset_query(set='C', fill='#77dd77'),
    upset_query(set='D', fill='#CD534CFF')
  ),
  base_annotations=list(
    'Intersection size'=intersection_size(
      text_mapping=aes(
        label=paste0(
          (round(
            !!get_size_mode('exclusive_intersection')/nrow(upset_data) * 100)),
          '%',
          '\n', '(',
          !!get_size_mode('exclusive_intersection'),
          ')'
        )
      ))
    + theme(
      # hide grid lines
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
    )
    + ylab('Combined Marker Positivity') 
    
  ),
  stripes=upset_stripes(
    geom=geom_segment(size=12),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the color
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(mapping=aes(y=..count../max(..count..)),
    )
    +scale_y_reverse(labels=scales::percent_format())
    +ylab("Marker Positivity")
    +geom_text(aes(label=..count..), hjust=0.5, stat="count", size=5)
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  ),
  sort_sets='descending',
  sort_intersections='descending'
)


# Scaled Venn diagram for binary readouts

venn_data <- upset_data

# Function to count the input when marker is positive
cts <- function(dataframe, condition) {
  df <- dataframe[condition,]
  n_rows <- nrow(df)
  return(n_rows)
}



library(eulerr)

v <- eulerr::euler(c(A= cts(venn_data, venn_data$A==1&venn_data$B==0&venn_data$C==0),
                     B=cts(venn_data, venn_data$A==0&venn_data$B==1&venn_data$C==0),
                     C=cts(venn_data, venn_data$A==0&venn_data$B=0&venn_data$C==1),
                     "A&B" =  cts(venn_data, venn_data$A==1&venn_data$B==1&venn_data$C==0),
                     "A&C" = cts(venn_data, venn_data$A ==1 & venn_data$B==0&venn_data$C==1),
                     "B&C" = cts(venn_data, venn_data$A==0&venn_data$B==1&venn_data$C==1),
                     "A&B&C" = cts(venn_data, venn_data$A==1&venn_data$B==1&venn_data$C==1),
                     "biomarker neg" = cts(venn_data, venn_data$A==0&venn_data$B==0&venn_data$C==0)),
                   shape = "ellipse") 


# Scaled venn plot should have same color codes as the upset
plot(v, key = TRUE, counts = TRUE, 
     quantities = list(type = c("counts", "percent"), 
                       font=1, round=5, cex=1),
     fill=c("#77dd77", "#CD534CFF", "#56B4E9"), 
     alpha = 0.3, edges=list(lty = 1), 
     factor_names = TRUE, labels=list(font=1, cex=1), 
     legend =TRUE) # With labels/markers


# Heatmap for continuous readouts/values

heatmap_data<- data

library(gplots)
library(ggplot2)
library(RColorBrewer)
library(viridis)

col<- viridis(100,alpha=1,begin = 0, end=1, option="D")
hclust.fun<- function(x) hclust(x, method="ward.D2")

# make sure all inputs are numeric
heatmap_data$A <- as.numeric(heatmap_data$A)
heatmap_data$B <- as.numeric(heatmap_data$B)
heatmap_data$C <- as.numeric(heatmap_data$C)
heatmap_data$D <- as.numeric(heatmap_data$D)


# create dendogram outside heatmap
mat<- as.matrix(heatmap_data)
row_clust <- hclust(dist(mat, method = 'euclidean'), method = 'ward.D2')
out <- heatmap.2( mat,col= col,
                  trace="none", density.info="none", scale="none",
                  main = "% positive cells-TICA", margins = c(7,10),
                  Rowv = as.dendrogram(row_clust))


# Summary Plot for continuous readouts/values

library(reshape2)
library(dplyr)

data_plot<- data
data_plot$CaseName <- rownames(data_plot) # CaseName is ID
data_plot<- melt(data_plot, id.vars=c("CaseName"), var=("MarkerName"))
colnames(data_plot) [3] <-  "memb(meanOD)_perc_pos"

# Creates a new column with the spacings (or location) or where to put the under bar within the ggplot object.

create_underbar <- function(data_plot){
  
  output <- data_plot %>%
    dplyr::mutate(underbar = dplyr::case_when(
      perc_positive_groups == 10 ~ -0.15,
      perc_positive_groups == 25 ~ -0.25,
      perc_positive_groups == 50 ~ -0.35,
      perc_positive_groups ==100 ~ -0.45
    ))
  
  return(output)
}

# Set the labels for the division of the 4 groups
labels <- c("<10%", "\u226510%", "\u226525%", "\u226550%")

# Set the numerical breaks within the `perc_positive_groups` column
breaks <- c(10,25,50,100)

# Select the colours for each group (HEX codes)
to_use <- c("#DC1310","#E69F00", "#F0E442", "#0072B2")

## add new column perc_positive_groups
data_plot$perc_positive_groups[data_plot$`memb(meanOD)_perc_pos` < 10] <- 10 
data_plot$perc_positive_groups[data_plot$`memb(meanOD)_perc_pos` >= 10 & data_plot$`memb(meanOD)_perc_pos` < 25] <- 25 
data_plot$perc_positive_groups[data_plot$`memb(meanOD)_perc_pos` >= 25 & data_plot$`memb(meanOD)_perc_pos` < 50] <- 50
data_plot$perc_positive_groups[data_plot$`memb(meanOD)_perc_pos` >= 50]  <- 100


# Select the columns required, drop all the NA values,
# and then convert the per_positive_groups column to the reversed order in breaks
plot.compact <- data_plot %>%
  dplyr::select(CaseName,MarkerName,perc_positive_groups) %>%
  tidyr::drop_na(perc_positive_groups) %>%
  dplyr::mutate(perc_positive_groups = factor(perc_positive_groups, levels = (breaks), ordered = TRUE))


# Number of samples per MarkerName
data_type_1 <- plot.compact %>% dplyr::count(MarkerName)

# Positions underneath the categorical heatmap where to but the coloured "underbars"
# See the create_underbar function below
data_type_2 <- plot.compact %>% dplyr::distinct(MarkerName, perc_positive_groups)
data_type_2 <- create_underbar(data = data_type_2)


data_type_3 <- plot.compact %>%
  dplyr::count(MarkerName, perc_positive_groups) %>%
  dplyr::group_by(MarkerName) %>%
  dplyr::mutate( total = sum(n)) 



b <- ggplot2::ggplot(plot.compact,
                     ggplot2::aes(x = MarkerName)) +
  ggplot2::geom_bar(ggplot2::aes(fill = perc_positive_groups),
                    position = ggplot2::position_fill(reverse = T)) +
  ggplot2::geom_text(ggplot2::aes(y = -0.05, label = paste0("n=", n)),
                     size = 10,
                     data = data_type_1) +
  ggplot2::geom_tile(ggplot2::aes(y = underbar ,width = 0.9, height = 0.1, fill = perc_positive_groups, color = NULL),
                     data = data_type_2) +
  ggplot2::geom_text(ggplot2::aes(y = ifelse(as.numeric(perc_positive_groups) == 1, -0.15, ifelse(as.numeric(perc_positive_groups) == 2, -0.25,
                                                                                                  ifelse(as.numeric(perc_positive_groups) == 3, -0.35,-0.45))),
                                  label = sprintf("%0.1f%%", n/total * 100)),
                     size = 10,
                     data = data_type_3) +
  ggplot2::scale_y_continuous(breaks = seq(0,1,0.25), labels = scales::label_percent()) +
  ggplot2::scale_fill_manual(values = to_use,labels = rev(labels), breaks = rev(breaks)) +
  ggplot2::labs(fill = "% positive cells",
                y = "fraction of samples",
                x = NULL,
                title = "xxx") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 legend.position = "top",
                 panel.spacing.x = ggplot2::unit(1, "cm")) +
  ggplot2::theme(text = ggplot2::element_text(size = 14))

b

