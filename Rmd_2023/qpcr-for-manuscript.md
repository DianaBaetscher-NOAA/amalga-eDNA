qPCR-analysis-for-manuscript
================
diana baetscher
2023-11-28

``` r
library(readxl)
library(patchwork)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggpubr)
```

## Clean analyses with the new qPCR data for manuscript

28 November 2023

``` r
# read in the extraction database metadata
metadata <- read_xlsx("../data/extraction_database.xlsx", sheet = "metadata")

# read in the combined qPCR results
#qpcr <- read_xlsx("../data/qPCR_results_combined.xlsx")

# organize color palette
mycols <- c("coral3", "cadetblue4")
```

### Metadata QC

``` r
# clean up the metadata format and variable names
frsh_meta <- metadata %>%
  mutate(tide = ifelse(`Sample Code` == "AM", "incoming", "outgoing")) %>%
  mutate(tide = ifelse(is.na(`Sample Code`) | `Sample Code` == "Creek", NA, tide)) %>%
  rename(sample = `Extraction ID`, location = `Sample Code`, label = `Sample Label`, type = `Sample_type`) %>%
  # distance needs to be a numeric variable without meters attached
  mutate(dist_grp = ifelse(location %in% c("PM", "AM"), label, NA))

# strip off the meters on dist_grp so we can order them properly
frsh_meta$dist_grp <- gsub("[^0-9.-]", "", frsh_meta$dist_grp)

# make another variable that is numeric distance (not character)
clean_meta <- frsh_meta %>%
  mutate(distance = as.numeric(dist_grp)) %>%
  # remove some of the unneeded variables
  dplyr::select(sample, location, label, type, tide, dist_grp, distance) %>%
  # add info to the creek samples
  mutate(label = ifelse(location == "Creek", "freshwater", label)) %>%
  mutate(label = ifelse(is.na(location), "blank", label))


clean_meta %>%
  write_csv("clean_meta.csv")
```

### qPCR data

## Combine data

Patterns across samples of inhibition?? Look at metadata.

Read in all combined data:

``` r
file_path <- "../data/qpcr/"

file_path %>%
  list.files() %>%
  .[str_detect(., ".xlsx")] -> xlsx_file_names

# Load everything into the Global Environment
xlsx_file_names %>%
  purrr::map(function(file_name){ # iterate through each file name
  
  read_xlsx(paste0(file_path, file_name), sheet = "Results", skip = 34)
  
}) -> df_list_read2 # Assign to a list

# combine all files into one df
# nothing was done to the quantities in this df
all_plates2021 <- bind_rows(df_list_read2, .id = "plate") %>%
    mutate(efficiency = -1 + 10^ (-1/Slope))

# calculate efficiency for all plates
all_plates2021 %>%
  group_by(`Target Name`) %>%
  summarise(E = mean(efficiency))
```

    ## # A tibble: 3 × 2
    ##   `Target Name`      E
    ##   <chr>          <dbl>
    ## 1 Moa_IPC        0.929
    ## 2 Oke_COI        0.907
    ## 3 <NA>          NA

``` r
# but I'm going to want to add a column that designates whether the sample is inhibited or not.
qpcr_data_w_metadata <- all_plates2021 %>%
  left_join(., clean_meta, by = c("Sample Name" = "sample"))

qpcr_data_w_metadata %>%
  filter(Task == "UNKNOWN" &
           CT != "Undetermined") %>%
  dplyr::select(`Sample Name`, `Ct Mean`, dist_grp, tide, distance) %>%
  unique()
```

    ## # A tibble: 276 × 5
    ##    `Sample Name` `Ct Mean` dist_grp tide     distance
    ##    <chr>             <dbl> <chr>    <chr>       <dbl>
    ##  1 e00056             28.4 160      incoming      160
    ##  2 e00056             25.6 160      incoming      160
    ##  3 e00128             28.6 2000     incoming     2000
    ##  4 e00190             28.7 1680     outgoing     1680
    ##  5 e00089             26.9 <NA>     <NA>           NA
    ##  6 e00125             27.4 1920     incoming     1920
    ##  7 e00189             27.2 1680     outgoing     1680
    ##  8 e00134             26.5 0        outgoing        0
    ##  9 e00134             25.5 0        outgoing        0
    ## 10 e00120             27.6 1840     incoming     1840
    ## # ℹ 266 more rows

Similar to the undiluted plates, I’ll want to read in all three results
documents and then analyze them in tandem.

### Diluted samples from 2021

``` r
diluted_file_path <- "../data/diluted_qpcr/"

diluted_file_path %>%
  list.files() %>%
  .[str_detect(., ".xlsx")] -> xlsx_file_names_diluted

# Load everything into the Global Environment
xlsx_file_names_diluted %>%
  purrr::map(function(file_name){ # iterate through each file name
  
  read_xlsx(paste0(diluted_file_path, file_name), sheet = "Results", skip = 34)
  
}) -> df_list_read3 # Assign to a list

diluted_plates <- bind_rows(df_list_read3, .id = "plate")

# multiple quantities by 10 because of 1:10 dilution factor
diluted_plates_10x <- diluted_plates %>%
  mutate(quant = Quantity*10) %>%
  mutate(inhibited = "yes")

# list of inhibited samples
inhibited_list <- diluted_plates_10x %>%
  filter(Task =="UNKNOWN") %>%
  dplyr::select(`Sample Name`) %>%
  unique()
```

### Combined data from inhibited samples and uninhibited samples

``` r
combined_data <- all_plates2021 %>%
  anti_join(., inhibited_list) %>%
  bind_rows(., diluted_plates_10x) %>% 
  filter(Task == "UNKNOWN" & 
           `Target Name` == "Oke_COI") %>%
  #select(`Sample Name`, `Quantity Mean`, quant) %>%
  mutate(quant = ifelse(is.na(quant), Quantity, quant)) %>%
  #filter(!is.na(quant)) %>%
  mutate(quant_per_ul = quant/2) %>%
  mutate(quant_per_L = quant/0.0016) # add a variable for the DNA copy number per L of seawater
```

    ## Joining with `by = join_by(`Sample Name`)`

``` r
# this maths out to 1L water x 0.4ml/5ml buffer x 2ul/100ul DNA = 1.6 mL

combined_data %>%
  write_csv("data_for_fig2.csv")

plotA_2021 <- combined_data %>%
  left_join(., clean_meta, by = c("Sample Name" = "sample")) %>%
  filter(!is.na(tide) &
           distance != 100) %>%
  ggplot(aes(x = distance, y = log10(quant_per_L), color = tide, shape = tide)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_smooth(method = "lm", se = F) +
  #facet_grid(rows = vars(tide)) +
  theme_bw() +
   scale_shape_manual(values = c(16, 18)) +
  scale_color_manual(values = mycols) +
  labs(x = "Distance from pens (m)",
       y = "DNA copies/L (log10)",
       color = "Tide",
       shape = "Tide") +
  theme(
    axis.title.x = element_text(margin = margin(t = 5)),
    #axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 5)),
        legend.position = c(0.8, 0.8)
    #legend.background = element_rect(color = "gray50")
  )

plotA_2021
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 198 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 198 rows containing missing values (`geom_point()`).

![](qpcr-for-manuscript_files/figure-gfm/qPCR-dataset-from-2021-1.png)<!-- -->

``` r
#ggsave("pdf_outputs/Amalga2021_qPCR_by_tide.png", width = 5, height = 4)
```

``` r
tmp2 <- combined_data %>%
  left_join(., clean_meta, by = c("Sample Name" = "sample"))%>%
  group_by(dist_grp,tide)%>%
  mutate(SampDist = as.numeric(as.factor(`Sample Name`)))

 tmp2 %>%
  filter(!is.na(tide) &
           distance != 100) %>%
  ggplot(aes(x = distance, y = log10(quant_per_L), color = tide)) +
  geom_point(aes(shape = as.factor(SampDist)),size = 2, alpha = 0.65) +
  geom_smooth(method = "lm", se = F) +
  #facet_grid(rows = vars(tide)) +
  theme_bw() +
  #scale_shape_manual(values = c(16, 18)) +
  scale_color_manual(values = mycols) +
  labs(x = "Distance from pens (m)",
       y = "DNA copies/L (log10)",
       color = "Tide",
       shape = "Water Sample") +
  theme(
    axis.title.x = element_text(margin = margin(t = 5)),
    #axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 5)),
    legend.position = c(0.8, 0.8)
    #legend.background = element_rect(color = "gray50")
  )
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 198 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 198 rows containing missing values (`geom_point()`).

![](qpcr-for-manuscript_files/figure-gfm/Pats-plot-for-water-replicates-1.png)<!-- -->

``` r
combined_data %>%
  left_join(., clean_meta, by = c("Sample Name" = "sample")) %>%
  filter(!is.na(tide) &
           distance != 100) %>%
  group_by(`Sample Name`, tide) %>%
  filter(distance < 10)
```

    ## # A tibble: 18 × 42
    ## # Groups:   Sample Name, tide [6]
    ##    plate  Well `Well Position` Omit  `Sample Name` `Target Name` Task   Reporter
    ##    <chr> <dbl> <chr>           <chr> <chr>         <chr>         <chr>  <chr>   
    ##  1 1        28 C4              false e00134        Oke_COI       UNKNO… FAM     
    ##  2 1        29 C5              false e00134        Oke_COI       UNKNO… FAM     
    ##  3 1        30 C6              false e00134        Oke_COI       UNKNO… FAM     
    ##  4 1        40 D4              false e00051        Oke_COI       UNKNO… FAM     
    ##  5 1        41 D5              false e00051        Oke_COI       UNKNO… FAM     
    ##  6 1        42 D6              false e00051        Oke_COI       UNKNO… FAM     
    ##  7 1        76 G4              false e00136        Oke_COI       UNKNO… FAM     
    ##  8 1        77 G5              false e00136        Oke_COI       UNKNO… FAM     
    ##  9 1        78 G6              false e00136        Oke_COI       UNKNO… FAM     
    ## 10 2        94 H10             false e00049        Oke_COI       UNKNO… FAM     
    ## 11 2        95 H11             false e00049        Oke_COI       UNKNO… FAM     
    ## 12 2        96 H12             false e00049        Oke_COI       UNKNO… FAM     
    ## 13 1        19 B7              false e00050        Oke_COI       UNKNO… FAM     
    ## 14 1        20 B8              false e00050        Oke_COI       UNKNO… FAM     
    ## 15 1        21 B9              false e00050        Oke_COI       UNKNO… FAM     
    ## 16 2        43 D7              false e00135        Oke_COI       UNKNO… FAM     
    ## 17 2        44 D8              false e00135        Oke_COI       UNKNO… FAM     
    ## 18 2        45 D9              false e00135        Oke_COI       UNKNO… FAM     
    ## # ℹ 34 more variables: Quencher <chr>, CT <chr>, `Ct Mean` <dbl>,
    ## #   `Ct SD` <dbl>, Quantity <dbl>, `Quantity Mean` <dbl>, `Quantity SD` <dbl>,
    ## #   `Automatic Ct Threshold` <chr>, `Ct Threshold` <dbl>,
    ## #   `Automatic Baseline` <chr>, `Baseline Start` <dbl>, `Baseline End` <dbl>,
    ## #   Comments <lgl>, `Y-Intercept` <dbl>, `R(superscript 2)` <dbl>, Slope <dbl>,
    ## #   `Amp Score` <dbl>, `Cq Conf` <dbl>, `Amp Status` <chr>, EXPFAIL <chr>,
    ## #   HIGHSD <chr>, NOAMP <chr>, efficiency <dbl>, AMPNC <chr>, quant <dbl>, …

At incoming vs. outgoing tides, the concentration of chum salmon DNA at
~1000 m looks different. Make a note of the LOD/LOQ in the figure
legend, but include detections below the LOQ on the plot.

Make a plot of the CV for each distance:

``` r
# plotB_2021 <- combined_data %>%
#   left_join(., clean_meta, by = c("Sample Name" = "sample")) %>%
#   filter(!is.na(tide) & 
#            !is.na(quant_per_ul) &
#            distance != 100) %>%
#   group_by(distance, tide) %>%
#   mutate(sd = sd(log10(quant_per_ul))) %>%
#   mutate(mean = mean(log10(quant_per_ul))) %>%
#   dplyr::select(distance, sd, mean, tide) %>%
#   unique() %>%
#   mutate(cv = sd/mean) %>%
#   ggplot(aes(x = distance, y = cv, color = tide, shape = tide)) +
#   geom_point(size = 2, alpha = 0.85) +
#   theme_bw() +
#    scale_shape_manual(values = c(16, 18)) +
#   scale_color_manual(values = mycols) +
#   labs(x = "Distance from pens (m)",
#        y = "CV",
#        color = "Tide",
#        shape = "Tide") +
#   theme(
#     axis.title.x = element_text(margin = margin(t = 5)),
#     axis.title.y = element_text(margin = margin(r = 5)),
#         #legend.position = c(0.8, 0.8)
#     legend.position = "none"
#   )

#ggsave("pdf_outputs/Amalga2021_qPCR_by_tideCV.png", width = 5, height = 4)
```

Which samples are inhibited?

``` r
fig2.a <- combined_data %>%
  left_join(., clean_meta, by = c("Sample Name" = "sample")) %>%
    filter(!is.na(tide) &
             distance != 100) %>%
  mutate(inhibited = ifelse(is.na(inhibited), "not inhibited", "inhibited")) %>%
  ggplot(aes(x = distance, y = log10(quant_per_L))) +
  geom_point(aes(color = inhibited), alpha = 0.8) +
  facet_grid(rows = vars(tide), labeller = label_both) +
  theme_bw() +
  labs(x = "Distance from pens (m)",
       y = "DNA copies/L (log10)") +
  theme(
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    legend.position = c(0.8, 0.9),
    #legend.background = element_rect(color = "gray50"),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm")
    #legend.spacing.y = unit(0.2, "cm")
  ) +
  scale_color_manual(values = c("darkblue", "darkgray"))

fig2.a
```

    ## Warning: Removed 198 rows containing missing values (`geom_point()`).

![](qpcr-for-manuscript_files/figure-gfm/inhibited-samples-plot-1.png)<!-- -->

``` r
#ggsave("pdf_outputs/Amalga2021_qPCR_inhibition_by_tide.png", width = 6, height = 5)
```

## 2022 data

``` r
fp <- "../data/qpcr/2022/"

fp %>%
  list.files() %>%
  .[str_detect(., ".xlsx")] -> xlsx_names

# Load everything into the Global Environment
xlsx_names %>%
  purrr::map(function(file_name){ # iterate through each file name
  
  read_xlsx(paste0(fp, file_name), sheet = "Results", skip = 34)
  
}) -> df_list1 # Assign to a list

plates2022 <- bind_rows(df_list1, .id = "plate")

# multiple quantities by 10 because of 1:10 dilution factor
# and 2 ul of diluted extract added to each well
plates2022_10x <- plates2022 %>%
  mutate(quant = (Quantity*10)/2) %>%
  mutate(quant_per_L = quant/0.0016) # add a variable for the DNA copy number per L of seawater
# this maths out to 1L water x 0.4ml/5ml buffer x 2ul/100ul DNA = 1.6 mL (because the 10x dilution was dealt with in the line above)

# list of samples
list2022 <- plates2022_10x %>%
  filter(Task =="UNKNOWN") %>%
  dplyr::select(`Sample Name`) %>%
  unique()

list2022
```

    ## # A tibble: 183 × 1
    ##    `Sample Name`
    ##    <chr>        
    ##  1 e01192       
    ##  2 e01256       
    ##  3 e01264       
    ##  4 e01193       
    ##  5 e01257       
    ##  6 e01265       
    ##  7 e01194       
    ##  8 e01258       
    ##  9 e01266       
    ## 10 e01195       
    ## # ℹ 173 more rows

``` r
# read in 2022 metadata 
meta2022 <- read_xlsx("../data/Amalga2022_cleanMetadata.xlsx")

meta2022
```

    ## # A tibble: 288 × 15
    ##    `Extraction ID` `Sample Code`   Project     `Sample Label`  transect     
    ##    <chr>           <chr>           <chr>       <chr>           <chr>        
    ##  1 e01192          0m_surface_AM   Amalga 2022 0m_surface_AM   perpendicular
    ##  2 e01193          0m_surface_AM   Amalga 2022 0m_surface_AM   perpendicular
    ##  3 e01194          0m_surface_AM   Amalga 2022 0m_surface_AM   perpendicular
    ##  4 e01195          0m_10m_depth_AM Amalga 2022 0m_10m_depth_AM perpendicular
    ##  5 e01196          0m_10m_depth_AM Amalga 2022 0m_10m_depth_AM perpendicular
    ##  6 e01197          0m_10m_depth_AM Amalga 2022 0m_10m_depth_AM perpendicular
    ##  7 e01198          0m_surface      Amalga 2022 0m_surface      perpendicular
    ##  8 e01199          0m_surface      Amalga 2022 0m_surface      perpendicular
    ##  9 e01200          0m_surface      Amalga 2022 0m_surface      perpendicular
    ## 10 e01201          500m_surface    Amalga 2022 500m_surface    perpendicular
    ## # ℹ 278 more rows
    ## # ℹ 10 more variables: distance_from_pens <dbl>, depth <dbl>,
    ## #   time_of_day <chr>, replicate <dbl>, Sample_Date <chr>,
    ## #   `Sampler Initials` <chr>, `Metadata Reference` <chr>,
    ## #   `Extraction Initials` <chr>, Ext_Date <dttm>, Ext_Plate <dttm>

``` r
qpcr_meta2022_combined <- plates2022_10x %>%
  left_join(., meta2022, by = c("Sample Name" = "Extraction ID")) %>%
  filter(Task == "UNKNOWN" &
          !is.na(distance_from_pens)) 

qpcr_meta2022_combined %>%
  filter(`Target Name` == "Oke_COI" &
           transect == "perpendicular") %>%
  ggplot(aes(x = distance_from_pens, y = log(quant))) +
  geom_point(aes(color = time_of_day), alpha = 0.5) +
  facet_grid(cols = vars(transect), rows = vars(depth), space = "free", labeller = label_both) +
  theme_bw()
```

    ## Warning: Removed 220 rows containing missing values (`geom_point()`).

![](qpcr-for-manuscript_files/figure-gfm/metadata-for-2022-samples-1.png)<!-- -->

``` r
list2022 %>%
  left_join(., meta2022, by = c("Sample Name" = "Extraction ID")) %>%
  filter(Sample_Date == "2022_05_05")
```

    ## # A tibble: 174 × 15
    ##    `Sample Name` `Sample Code`     Project     `Sample Label`    transect     
    ##    <chr>         <chr>             <chr>       <chr>             <chr>        
    ##  1 e01192        0m_surface_AM     Amalga 2022 0m_surface_AM     perpendicular
    ##  2 e01264        500m_10m_depth_AM Amalga 2022 500m_10m_depth_AM perpendicular
    ##  3 e01193        0m_surface_AM     Amalga 2022 0m_surface_AM     perpendicular
    ##  4 e01265        500m_10m_depth_AM Amalga 2022 500m_10m_depth_AM perpendicular
    ##  5 e01194        0m_surface_AM     Amalga 2022 0m_surface_AM     perpendicular
    ##  6 e01258        500m_surface_AM   Amalga 2022 500m_surface_AM   perpendicular
    ##  7 e01266        500m_10m_depth_AM Amalga 2022 500m_10m_depth_AM perpendicular
    ##  8 e01195        0m_10m_depth_AM   Amalga 2022 0m_10m_depth_AM   perpendicular
    ##  9 e01259        500m_surface_AM   Amalga 2022 500m_surface_AM   perpendicular
    ## 10 e01267        1000m_surface_AM  Amalga 2022 1000m_surface_AM  perpendicular
    ## # ℹ 164 more rows
    ## # ℹ 10 more variables: distance_from_pens <dbl>, depth <dbl>,
    ## #   time_of_day <chr>, replicate <dbl>, Sample_Date <chr>,
    ## #   `Sampler Initials` <chr>, `Metadata Reference` <chr>,
    ## #   `Extraction Initials` <chr>, Ext_Date <dttm>, Ext_Plate <dttm>

``` r
# 2021 data - including the inhibited samples rerun as dilutions
df2021 <- combined_data %>%
  filter(`Target Name` == "Oke_COI" &
           Task == "UNKNOWN") %>%
  dplyr::select(`Sample Name`, CT, quant_per_ul, quant_per_L)
  
  
# 2022 data
df2022 <- plates2022_10x %>%
  filter(`Target Name` == "Oke_COI" &
           Task == "UNKNOWN") %>%
  rename(quant_per_ul = quant) %>%
  dplyr::select(`Sample Name`, CT, quant_per_ul, quant_per_L)


# combine those
combo_df <- bind_rows(df2021, df2022) %>%
  rename(sample = `Sample Name`)
```

``` r
slim_m2021 <- clean_meta %>%
  mutate(time_of_day = location) %>%
  mutate(depth = 0) %>%
  mutate(transect = "perpendicular") %>%
  rename(distance_from_pens = distance) %>%
  dplyr::select(-location, -type, -dist_grp) %>%
  mutate(sample_date = "2021_05_10")

slim_m2022 <- meta2022 %>%
  rename(sample = `Extraction ID`, label = `Sample Code`, sample_date = Sample_Date) %>%
  dplyr::select(sample, label, distance_from_pens, time_of_day, depth, transect, sample_date) %>%
  mutate(tide = ifelse(time_of_day == "AM", "outgoing", "incoming")) 
  

# combine metadata from both years
metadata_both_years <- bind_rows(slim_m2021, slim_m2022) 
```

``` r
merged_df_both_years <- combo_df %>%
  left_join(., metadata_both_years) %>%
  filter(sample_date %in% c("2021_05_10", "2022_05_05")) %>%
  mutate(year = ifelse(sample_date == "2021_05_10", "2021", "2022")) %>%
  filter(distance_from_pens != 100) %>%
    mutate(log_dna_L = log10(quant_per_L)) #### Here, I've changed the per ul (rxn volume) to per L of seawater sampled.
```

    ## Joining with `by = join_by(sample)`

``` r
merged_df_both_years %>%
  filter(year == "2022" & depth == 0)
```

    ## # A tibble: 162 × 13
    ##    sample CT           quant_per_ul quant_per_L label   tide  distance_from_pens
    ##    <chr>  <chr>               <dbl>       <dbl> <chr>   <chr>              <dbl>
    ##  1 e01192 32.188362            278.     173888. 0m_sur… outg…                  0
    ##  2 e01192 31.80571             356.     222709. 0m_sur… outg…                  0
    ##  3 e01192 32.261612            265.     165843. 0m_sur… outg…                  0
    ##  4 e01193 33.708157            104.      65078. 0m_sur… outg…                  0
    ##  5 e01193 32.88484             177.     110831. 0m_sur… outg…                  0
    ##  6 e01193 32.897076            176.     109958. 0m_sur… outg…                  0
    ##  7 e01194 30.890171            644.     402594. 0m_sur… outg…                  0
    ##  8 e01194 30.891844            643.     402159. 0m_sur… outg…                  0
    ##  9 e01194 30.88471             646.     404019. 0m_sur… outg…                  0
    ## 10 e01258 Undetermined          NA          NA  500m_s… outg…                500
    ## # ℹ 152 more rows
    ## # ℹ 6 more variables: time_of_day <chr>, depth <dbl>, transect <chr>,
    ## #   sample_date <chr>, year <chr>, log_dna_L <dbl>

``` r
merged_df_both_years %>%
  filter(year == "2022") %>%
  group_by(sample, tide, depth) %>%
  tally() %>%
  group_by(tide, depth) %>%
  tally()
```

    ## # A tibble: 6 × 3
    ## # Groups:   tide [2]
    ##   tide     depth     n
    ##   <chr>    <dbl> <int>
    ## 1 incoming     0    27
    ## 2 incoming     5    27
    ## 3 incoming    10    27
    ## 4 outgoing     0    27
    ## 5 outgoing     5    27
    ## 6 outgoing    10    27

``` r
merged_df_both_years %>%
  filter(year == 2021) %>%
  group_by(sample, tide) %>%
  tally() %>%
  group_by(tide) %>%
  tally()
```

    ## # A tibble: 2 × 2
    ##   tide         n
    ##   <chr>    <int>
    ## 1 incoming    78
    ## 2 outgoing    78

Okay, that’s kind of useful. Very little to talk about with the 2022
samples. There’s one outlier for the incoming tide at 0m that I don’t
like - it would be nice to take a look at the replicates.

``` r
merged_df_both_years %>%
  filter(sample_date == "2022_05_05") %>%
  mutate(detection = ifelse(CT == "Undetermined", "no", "yes")) %>%
  group_by(detection) %>%
  tally()
```

    ## # A tibble: 2 × 2
    ##   detection     n
    ##   <chr>     <int>
    ## 1 no          425
    ## 2 yes          61

So many non-detections!

``` r
merged_df_both_years %>%
  filter(sample_date == "2021_05_10") %>%
  mutate(detection = ifelse(CT == "Undetermined", "no", "yes")) %>%
  group_by(detection) %>%
  tally()
```

    ## # A tibble: 2 × 2
    ##   detection     n
    ##   <chr>     <int>
    ## 1 no          198
    ## 2 yes         269

## Stats for distance & tide

Are data normally distributed?

``` r
for_mod2021 <- merged_df_both_years %>%
  filter(sample_date == "2021_05_10" &
           distance_from_pens != 100) # bec this is between the creek and the pens




ggdensity(for_mod2021$log_dna_L)
```

    ## Warning: Removed 198 rows containing non-finite values (`stat_density()`).

![](qpcr-for-manuscript_files/figure-gfm/test-2021-data-for-normality-1.png)<!-- -->

``` r
ggqqplot(for_mod2021$log_dna_L)
```

    ## Warning: Removed 198 rows containing non-finite values (`stat_qq()`).

    ## Warning: Removed 198 rows containing non-finite values (`stat_qq_line()`).
    ## Removed 198 rows containing non-finite values (`stat_qq_line()`).

![](qpcr-for-manuscript_files/figure-gfm/test-2021-data-for-normality-2.png)<!-- -->

Those visual inspection methods make the data appear normally
distributed, but it’s good to test this explicitly.

``` r
shapiro.test(for_mod2021$log_dna_L)
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  for_mod2021$log_dna_L
    ## W = 0.98787, p-value = 0.02323

Based on this, the data are not significantly different from normal.

``` r
# model 1: copy number and distance
mod1 <- lm(data = for_mod2021, log_dna_L ~ distance_from_pens)

# model 2: distance and tide
mod2 <- lm(data = for_mod2021, log_dna_L ~ distance_from_pens + tide)

# now compare the models with anova
anova(mod1, mod2)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: log_dna_L ~ distance_from_pens
    ## Model 2: log_dna_L ~ distance_from_pens + tide
    ##   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
    ## 1    267 79.205                                  
    ## 2    266 56.685  1     22.52 105.68 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(mod1)
```

    ## 
    ## Call:
    ## lm(formula = log_dna_L ~ distance_from_pens, data = for_mod2021)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.41502 -0.34342  0.00262  0.41874  1.34364 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         6.003e+00  6.014e-02   99.81   <2e-16 ***
    ## distance_from_pens -1.784e-03  7.731e-05  -23.08   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5447 on 267 degrees of freedom
    ##   (198 observations deleted due to missingness)
    ## Multiple R-squared:  0.6661, Adjusted R-squared:  0.6649 
    ## F-statistic: 532.7 on 1 and 267 DF,  p-value: < 2.2e-16

``` r
summary(mod2)
```

    ## 
    ## Call:
    ## lm(formula = log_dna_L ~ distance_from_pens + tide, data = for_mod2021)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.15499 -0.29317  0.03599  0.32653  1.68978 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         6.2893739  0.0581049  108.24   <2e-16 ***
    ## distance_from_pens -0.0018164  0.0000656  -27.69   <2e-16 ***
    ## tideoutgoing       -0.5814644  0.0565630  -10.28   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4616 on 266 degrees of freedom
    ##   (198 observations deleted due to missingness)
    ## Multiple R-squared:  0.7611, Adjusted R-squared:  0.7593 
    ## F-statistic: 423.6 on 2 and 266 DF,  p-value: < 2.2e-16

``` r
library(AICcmodavg)

# list of models
models <- list(mod1, mod2)

mod.names <- c("dist", "dist.tide")

# calculate AIC of each model
AICcmodavg::aictab(cand.set = models, modnames = mod.names)
```

    ## 
    ## Model selection based on AICc:
    ## 
    ##           K   AICc Delta_AICc AICcWt Cum.Wt      LL
    ## dist.tide 4 352.65       0.00      1      1 -172.25
    ## dist      3 440.58      87.93      0      1 -217.24

``` r
merged_df_both_years %>%
  filter(year == "2022" &
           distance_from_pens == 0) %>%
  group_by(tide) %>%
  filter(!is.na(quant_per_L)) %>%
  summarise(mean(quant_per_L))
```

    ## # A tibble: 2 × 2
    ##   tide     `mean(quant_per_L)`
    ##   <chr>                  <dbl>
    ## 1 incoming              47266.
    ## 2 outgoing             161525.

``` r
fig2.d <- merged_df_both_years %>%
 filter(!is.na(distance_from_pens) &
           year == 2021 &
          distance_from_pens != 100) %>%
  group_by(tide, distance_from_pens) %>%
  mutate(detection = ifelse(is.na(quant_per_L), 0, 1)) %>%
  summarise(n_detections = sum(detection)/9) %>%
  ggplot(aes(x = distance_from_pens, y = n_detections, fill = tide)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(tide), labeller = label_both) +
  theme_bw() +
  scale_shape_manual(values = c(1,2)) +
  scale_fill_manual(values = mycols) +
  labs(x = "Distance from pens (m)",
       y = "Proportion of detections") +
  theme(
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    legend.position = "none"
  ) 
```

    ## `summarise()` has grouped output by 'tide'. You can override using the
    ## `.groups` argument.

``` r
#ggsave("pdf_outputs/eDNA_detections2021.png", width = 6, height = 4)
```

``` r
# combine the pieces from the other plots
(fig2.a | fig2.d)/plotA_2021 + plot_annotation(tag_levels = "A")
```

    ## Warning: Removed 198 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 198 rows containing non-finite values (`stat_smooth()`).
    ## Removed 198 rows containing missing values (`geom_point()`).

![](qpcr-for-manuscript_files/figure-gfm/make-figure-2-1.png)<!-- -->

``` r
ggsave("pdf_outputs/Figure2_tri.png", width = 9, height = 8)
```

    ## Warning: Removed 198 rows containing missing values (`geom_point()`).

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 198 rows containing non-finite values (`stat_smooth()`).
    ## Removed 198 rows containing missing values (`geom_point()`).

#### Breakpoint analysis for detections in 2021

Based on review feedback, I’ll analyze the proportion of detection data
from 2021 using a breakpoint (or segmented regression) analysis.

following the analysis here: <https://rpubs.com/MarkusLoew/12164>

``` r
detections2021 <- merged_df_both_years %>%
 filter(!is.na(distance_from_pens) &
           year == 2021 &
          distance_from_pens != 100) %>%
  group_by(tide, distance_from_pens) %>%
  mutate(detection = ifelse(is.na(quant_per_L), 0, 1)) %>%
  summarise(n_detections = sum(detection)/9) 
```

    ## `summarise()` has grouped output by 'tide'. You can override using the
    ## `.groups` argument.

``` r
# plot the detections
p <- ggplot(detections2021, aes(x = distance_from_pens, y = n_detections, color = tide)) + geom_line()

p
```

![](qpcr-for-manuscript_files/figure-gfm/detections-for-breakpoint-initial-setup-1.png)<!-- -->

``` r
# create a linear model
detect.lm <- lm(n_detections ~ distance_from_pens + tide, data = detections2021)
summary(detect.lm)
```

    ## 
    ## Call:
    ## lm(formula = n_detections ~ distance_from_pens + tide, data = detections2021)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.5183 -0.1384 -0.0384  0.1458  0.4883 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.185e+00  6.939e-02  17.075  < 2e-16 ***
    ## distance_from_pens -5.610e-04  5.291e-05 -10.603 2.78e-14 ***
    ## tideoutgoing       -9.829e-02  6.350e-02  -1.548    0.128    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2289 on 49 degrees of freedom
    ## Multiple R-squared:  0.7009, Adjusted R-squared:  0.6887 
    ## F-statistic: 57.41 on 2 and 49 DF,  p-value: 1.44e-13

``` r
# a linear model with data for the part after 1000 m
detect.lm2 <- lm(n_detections ~ distance_from_pens + tide, data = detections2021[detections2021$distance_from_pens > 1000, ])
summary(detect.lm2)
```

    ## 
    ## Call:
    ## lm(formula = n_detections ~ distance_from_pens + tide, data = detections2021[detections2021$distance_from_pens > 
    ##     1000, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.46398 -0.11935  0.01862  0.13172  0.38828 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)         1.7197802  0.2328204   7.387 1.64e-07 ***
    ## distance_from_pens -0.0009234  0.0001476  -6.258 2.19e-06 ***
    ## tideoutgoing       -0.0769231  0.0883365  -0.871    0.393    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2252 on 23 degrees of freedom
    ## Multiple R-squared:  0.6345, Adjusted R-squared:  0.6027 
    ## F-statistic: 19.96 on 2 and 23 DF,  p-value: 9.416e-06

``` r
# Extract te coefficients from the overall model
my.coef <- coef(detect.lm)

# add the regression line to the graph
# setting the aesthetics to a constant - this provides a name that we can reference later when we add additional layers
p <- p + geom_abline(intercept = my.coef[1], 
                         slope = my.coef[2], 
                     aes(colour = "overall"))
```

    ## Warning: `geom_abline()`: Ignoring `mapping` because `slope` and/or `intercept` were
    ## provided.

``` r
p
```

![](qpcr-for-manuscript_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# -------------------
# analyse breakpoints
# -------------------
# http://cran.r-project.org/doc/Rnews/Rnews_2008-1.pdf
library(segmented)
```

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## The following object is masked from 'package:patchwork':
    ## 
    ##     area

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

``` r
# have to provide estimates for breakpoints.
# after looking a the data, 
my.seg <- segmented(detect.lm, 
                    seg.Z = ~ distance_from_pens, 
                    psi = list(distance_from_pens = c(1200)))

# When not providing estimates for the breakpoints "psi = NA" can be used.
# The number of breakpoints that will show up is not defined
#my.seg <- segmented(my.lm, 
#                    seg.Z = ~ DistanceMeters, 
#                    psi = NA)

# display the summary
summary(my.seg)
```

    ## 
    ##  ***Regression Model with Segmented Relationship(s)***
    ## 
    ## Call: 
    ## segmented.lm(obj = detect.lm, seg.Z = ~distance_from_pens, psi = list(distance_from_pens = c(1200)))
    ## 
    ## Estimated Break-Point(s):
    ##                          Est.  St.Err
    ## psi1.distance_from_pens  960 140.156
    ## 
    ## Meaningful coefficients of the linear terms:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)            0.9890137  0.0783636  12.621   <2e-16 ***
    ## distance_from_pens    -0.0001447  0.0001297  -1.115    0.271    
    ## tideoutgoing          -0.0982906  0.0549232  -1.790    0.080 .  
    ## U1.distance_from_pens -0.0007872  0.0001835  -4.290       NA    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.198 on 47 degrees of freedom
    ## Multiple R-Squared: 0.7853,  Adjusted R-squared: 0.7671 
    ## 
    ## Boot restarting based on 6 samples. Last fit:
    ## Convergence attained in 2 iterations (rel. change 2.6511e-08)

``` r
merged_df_both_years %>%
  filter(transect == "perpendicular" &
           depth == 0) %>%
  ggplot(aes(x = distance_from_pens, y = log_dna_L, color = year)) +
  geom_point() +
  facet_grid(rows = vars(tide)) +
  theme_bw() +
  scale_color_manual(values = c("darksalmon", "dodgerblue")) +
  labs(x = "Distance from pens (m)",
        y = "DNA copies/L (log10)",
       color = "Sampling year") +
  theme(
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5))
  ) 
```

    ## Warning: Removed 256 rows containing missing values (`geom_point()`).

![](qpcr-for-manuscript_files/figure-gfm/make-si-figure-comparing-conc-across-years-1.png)<!-- -->

``` r
ggsave("pdf_outputs/SI_figure_S5_conc_across_years.png", width = 6, height = 4)
```

    ## Warning: Removed 256 rows containing missing values (`geom_point()`).

``` r
df_for_3d <- merged_df_both_years %>%
  filter(year == "2022") %>%
  mutate(distance_from_pens = ifelse(str_detect(transect, "_B"), -distance_from_pens, distance_from_pens)) %>%
  mutate(x_dist = ifelse(transect == "perpendicular", distance_from_pens, NA)) %>%
  mutate(y_dist = ifelse(str_detect(transect,"parallel"), distance_from_pens, NA)) %>%
  mutate(y_dist = ifelse(is.na(y_dist), 0, y_dist)) %>%
  mutate(x_dist = ifelse(is.na(x_dist), 1000, x_dist)) %>%
  mutate(log_quant = log10(quant_per_ul)) %>%
  mutate(depth = ifelse(depth != 0, -1*depth, depth))
```

``` r
# LOD vs. LOQ
df_for_3d %>%
  filter(quant_per_ul > 38) %>%
  group_by(tide, distance_from_pens) %>%
  summarise(mean(quant_per_ul))
```

    ## `summarise()` has grouped output by 'tide'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 3 × 3
    ## # Groups:   tide [2]
    ##   tide     distance_from_pens `mean(quant_per_ul)`
    ##   <chr>                 <dbl>                <dbl>
    ## 1 incoming                  0                129. 
    ## 2 incoming                500                 44.6
    ## 3 outgoing                  0                366.

``` r
# I want to explicitly plot NAs on here to show where we sampled but didn't have data
# df_for_3d %>%
#   mutate(depth = -1*depth) %>%
#   filter(!is.na(depth)) %>%
#   ggplot(aes(x = x_dist, y = y_dist, size = log_quant, color = tide, shape = tide)) +
#   geom_jitter(alpha = 0.5, width = 100, height = 100) +
#   theme_bw() +
#   facet_grid(rows = vars(depth), labeller = label_both) +
#   labs(x = "Perpendicular transect (meters from pens)",
#        y = "Parallel transect (SE-NW, meters from 1000 m midpoint)",
#        color = "Tide",
#        shape = "Tide",
#        size = "DNA copies/uL\n(log10)") +
#   theme(
#     axis.title.x = element_text(margin = margin(t = 10)),
#     axis.title.y = element_text(margin = margin(r = 10))
#   ) +
#   scale_shape_manual(values = c(16, 18), na.value = 1) +
#   scale_color_manual(values = mycols) +
#   guides(color = guide_legend(override.aes = list(size = 5)))
# 
# ggsave("pdf_outputs/depth_transects2022.png", width = 6, height = 7)
```

``` r
library(viridis)
```

    ## Loading required package: viridisLite

``` r
# I want to explicitly plot NAs on here to show where we sampled but didn't have data
tmp <- df_for_3d %>%
  mutate(depth = -1*depth) %>%
  filter(!is.na(depth))

tmp %>% write_csv("amalga_tmp_df.csv")

## modifications from Pat for manual jittering
tmp<-tmp %>%
       mutate(yJit = rep(c(-100,0,100),nrow(tmp)/3))%>%
       group_by(label)%>%
       mutate(xJit = as.numeric(as.factor(label))+as.numeric(as.factor(sample))-2)%>%
       mutate(xJit = recode(xJit,
                            `1` = -100,
                            `2` = 100))
tmp %>%
  ggplot() +
  geom_point(aes(x = x_dist+xJit, y = y_dist+yJit, color = log_dna_L), size = 1, alpha = 0.9, width = 110, height = 90) +
  geom_point(data = tmp[is.na(tmp$log_dna_L),], aes(x_dist, y_dist), color = "gray", shape = 1, size = 3) +
  theme_bw() +
  facet_grid(rows = vars(depth), cols = vars(tide), labeller = label_both) +
  scale_color_viridis_c(option="plasma", na.value = NA, direction = -1) +
  labs(x = "Meters from pens along perpendicular transect",
       y = "Meters along parallel transect from midpoint (SE-NW)",
       color = "DNA copies/L\n(log10)") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )
```

    ## Warning in geom_point(aes(x = x_dist + xJit, y = y_dist + yJit, color =
    ## log_dna_L), : Ignoring unknown parameters: `width` and `height`

    ## Warning: Removed 425 rows containing missing values (`geom_point()`).

![](qpcr-for-manuscript_files/figure-gfm/explicitly-plot-NAs-1.png)<!-- -->

``` r
ggsave("pdf_outputs/depth_transects2022_v2.png", width = 6, height = 5)
```

    ## Warning: Removed 425 rows containing missing values (`geom_point()`).

Quick look at the dispersion across the parallel transect:

``` r
df_for_3d %>%
  filter(#depth == 0 &
           x_dist == 1000 #&
           #!is.na(log_quant)
             ) %>%
  group_by(depth, tide) %>%
  tally()
```

    ## # A tibble: 6 × 3
    ## # Groups:   depth [3]
    ##   depth tide         n
    ##   <dbl> <chr>    <int>
    ## 1   -10 incoming    45
    ## 2   -10 outgoing    45
    ## 3    -5 incoming    45
    ## 4    -5 outgoing    45
    ## 5     0 incoming    45
    ## 6     0 outgoing    45

``` r
  # ggplot(aes(x = y_dist, y = log_dna_ul, color = tide)) +
  # geom_point(size = 2, alpha = 0.5)
```

Relationship between qPCR detections and depth

``` r
# df_for_3d %>%
#   mutate(depth = -1*depth) %>%
#   filter(!is.na(depth) &
#            transect == "perpendicular") %>%
#   ggplot(aes(x = depth, y = log_dna_L, color = log_quant)) +
#   geom_point() +
#   facet_grid(rows = vars(tide), cols = vars(transect)) +
#   scale_x_continuous(breaks = c(0, 5, 10)) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(margin = margin(t = 10)),
#     axis.title.y = element_text(margin = margin(r = 10))
#   ) +
#     scale_color_gradient(na.value = "gray") +
#   guides(color = guide_legend(override.aes = list(size = 5))) +
#    labs(x = "Depth (m)",
#        y = "DNA copies/uL\n(log10)",
#        color = "Tide",
#        shape = "Tide",
#        size = "DNA copies/L\n(log10)")
```

Is it worth fitting a linear model to this?

``` r
# dis22p <- df_for_3d %>%
#   #mutate(depth = -1*depth) %>%
#   filter(!is.na(depth) &
#            transect == "perpendicular") %>%
#   ggplot(aes(x = x_dist, y = depth, size = log_quant, color = tide, shape = tide)) +
#   geom_jitter(alpha = 0.5, width = 70, height = 0.5) +
#   facet_grid(rows = vars(tide), cols = vars(transect)) +
#   #geom_smooth(method = "lm", se = F) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(margin = margin(t = 10)),
#     axis.title.y = element_text(margin = margin(r = 10))
#   ) +
#   scale_shape_manual(values = c(16, 18)) +
#   scale_color_manual(values = mycols) +
#   scale_y_continuous(breaks = c(0, -5, -10)) +
#    guides(color = guide_legend(override.aes = list(size = 5))) +
#    labs(y = "Depth (m)",
#        x = "Distance from pens (m)",
#        color = "Tide",
#        shape = "Tide",
#        size = "DNA copies/uL\n(log10)")
# 
# ggsave("pdf_outputs/depth_by_distance.png", width = 6, height = 5)
```

``` r
# df_for_3d %>%
#   #mutate(depth = -1*depth) %>%
#   mutate(transect = ifelse(distance_from_pens == 1000 & transect == "perpendicular", "parallel", transect)) %>%
#   mutate(y_dist = ifelse(distance_from_pens == 1000 & transect == "parallel", 0, y_dist)) %>%
#   filter(!is.na(depth) &
#            str_detect(transect, "parallel")) %>%
#   ggplot(aes(x = y_dist, y = depth, size = log_quant, color = tide, shape = tide)) +
#   geom_jitter(alpha = 0.5, width = 25, height = 0.3) +
#   #facet_grid(rows = vars(tide), cols = vars(transect)) +
#   #geom_smooth(method = "lm", se = F) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_text(margin = margin(t = 10)),
#     axis.title.y = element_text(margin = margin(r = 10))
#   ) +
#   scale_shape_manual(values = c(16, 18)) +
#   scale_color_manual(values = mycols) +
#   scale_y_continuous(breaks = c(0, -5, -10)) +
#    guides(color = guide_legend(override.aes = list(size = 5))) +
#    labs(y = "Depth (m)",
#        x = "Distance from pens (m)",
#        color = "Tide",
#        shape = "Tide",
#        size = "DNA copies/uL\n(log10)")
# 
# ggsave("pdf_outputs/depth_by_distance_parallel.png", width = 6, height = 5)
```

% of detections in surface samples vs. at depth?

### Depth Stats - ANOVAs

1.  Test whether there is statistically more DNA at the surface compared
    to the 5m and 10m depths
2.  and then test whether incoming/outgoing tide makes a difference

``` r
for_anova <- merged_df_both_years %>%
  filter(year == 2022 &
           !is.na(depth) &
           !is.na(quant_per_ul) &
           distance_from_pens == 0) %>%
  group_by(depth, tide) 

for_anova %>%
  summarise(sd = sd(log_dna_L))
```

    ## `summarise()` has grouped output by 'depth'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 6 × 3
    ## # Groups:   depth [3]
    ##   depth tide        sd
    ##   <dbl> <chr>    <dbl>
    ## 1     0 incoming 0.302
    ## 2     0 outgoing 0.287
    ## 3     5 incoming 0.242
    ## 4     5 outgoing 0.118
    ## 5    10 incoming 1.41 
    ## 6    10 outgoing 0.378

``` r
for_anova %>%
  filter(quant_per_ul > 38)
```

    ## # A tibble: 17 × 13
    ## # Groups:   depth, tide [2]
    ##    sample CT        quant_per_ul quant_per_L label      tide  distance_from_pens
    ##    <chr>  <chr>            <dbl>       <dbl> <chr>      <chr>              <dbl>
    ##  1 e01192 32.188362        278.      173888. 0m_surfac… outg…                  0
    ##  2 e01192 31.80571         356.      222709. 0m_surfac… outg…                  0
    ##  3 e01192 32.261612        265.      165843. 0m_surfac… outg…                  0
    ##  4 e01193 33.708157        104.       65078. 0m_surfac… outg…                  0
    ##  5 e01193 32.88484         177.      110831. 0m_surfac… outg…                  0
    ##  6 e01193 32.897076        176.      109958. 0m_surfac… outg…                  0
    ##  7 e01194 30.890171        644.      402594. 0m_surfac… outg…                  0
    ##  8 e01194 30.891844        643.      402159. 0m_surfac… outg…                  0
    ##  9 e01194 30.88471         646.      404019. 0m_surfac… outg…                  0
    ## 10 e01312 33.755558         98.1      61310. 0m_surfac… inco…                  0
    ## 11 e01312 34.813236         49.8      31117. 0m_surfac… inco…                  0
    ## 12 e01313 33.640163        106.       66018. 0m_surfac… inco…                  0
    ## 13 e01313 33.24714         136.       84939. 0m_surfac… inco…                  0
    ## 14 e01313 34.739902         52.2      32615. 0m_surfac… inco…                  0
    ## 15 e01310 32.73228         189.      118162. 0m_surfac… inco…                  0
    ## 16 e01310 32.447186        227.      141862. 0m_surfac… inco…                  0
    ## 17 e01310 32.888504        171.      106899. 0m_surfac… inco…                  0
    ## # ℹ 6 more variables: time_of_day <chr>, depth <dbl>, transect <chr>,
    ## #   sample_date <chr>, year <chr>, log_dna_L <dbl>

171 samples per depth from 2022, but the number of detections per depth
is much smaller.

multiple comparisons: every depth x tide compared to every other depth x
tide

``` r
# set-up df
for_anova$depth <- as.factor(for_anova$depth)
for_anova$tide <- as.factor(for_anova$tide)

# form
# aov(response ~ factor A/ factor B)

depth.tide <- aov(for_anova$log_dna_L ~ for_anova$tide / for_anova$depth, data = for_anova)
summary(depth.tide)
```

    ##                                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## for_anova$tide                  1  2.137  2.1371   12.70  0.00174 ** 
    ## for_anova$tide:for_anova$depth  4 11.062  2.7655   16.43 2.28e-06 ***
    ## Residuals                      22  3.702  0.1683                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
tukey.depth.tide <- TukeyHSD(depth.tide, conf.level = 0.95)
as.data.frame(tukey.depth.tide$`for_anova$tide:for_anova$depth`) %>%
  arrange(`p adj`) %>%
  rownames_to_column(var = "tide.depth") %>%
  write_csv("tukey_depth_tide.csv")
```

``` r
plot(TukeyHSD(depth.tide, conf.level = 0.95), las = 2)
```

![](qpcr-for-manuscript_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](qpcr-for-manuscript_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
box_plot_df <- for_anova %>%
  mutate(log_dna_L = ifelse(log_dna_L < 0, 0, log_dna_L))
  

d.boxp <- ggplot(box_plot_df, aes(x=factor(depth), y=log_dna_L, fill=tide)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("coral3", "cadetblue4")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
    #legend.position = c(0.85,0.85)
  ) +
   labs(x = "Depth (m)",
       y = "DNA copies/L(log10)",
       fill = "Tide")

ggsave("pdf_outputs/depth_tide_boxplot.png", height = 3, width = 5)
```

``` r
# make figure 4
#dis22p + d.boxp + plot_annotation(tag_levels = "A") + plot_layout(ncol = 1, heights = c(3,2))

#ggsave("pdf_outputs/depth_tide_anova_plot.png", width = 6, height = 7)
```

### Depth eDNA below the LOQ

Because nearly all non-surface samples were below the LOQ, let’s switch
this to an analysis using detections rather than DNA concentration.

``` r
for_anova %>%
  mutate(detection = ifelse(quant_per_ul > 0, 1,0)) %>%
  group_by(depth, tide) %>%
  summarise(n_detections = sum(detection)) %>%
  mutate(prop_detections = n_detections/9) %>%
  ggplot(aes(x = depth, y = prop_detections, fill = tide)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(width=0.6)) +
  theme_minimal() +
  scale_fill_manual(values = c("coral3", "cadetblue4")) +
  theme_bw() +
  theme(
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = c(0.75,0.75)) +
   labs(x = "Depth (m)",
       y = "Proportion of eDNA detections",
       fill = "Tide") +
   scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0))
```

    ## `summarise()` has grouped output by 'depth'. You can override using the
    ## `.groups` argument.

![](qpcr-for-manuscript_files/figure-gfm/detections-by-tide-at-any-quantity-1.png)<!-- -->

``` r
ggsave("pdf_outputs/depth_detections_by_tide.png", width = 3.2, height = 2.8)
```
