# This is an R script
# Lines starting with a hash are comments. They are not executed as 
# commands

# To run a command, put your cursor on the line you want to run
# and do CTRL-ENTER. You can also use the 'Run' button to run a command.
# The command is executed in the console (below) and any output 
# will appear there
# the > in the console is the prompt. 


# load packages needed
library("phyloseq")
library("tidyverse")
library("ggvenn")
# The output from these commands will be:
### ── Attaching packages ──────────────────────────────────────── tidyverse 1.3.2 ──
### ✔ ggplot2 3.3.6      ✔ purrr   0.3.4 
### ✔ tibble  3.1.8      ✔ dplyr   1.0.10
### ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
### ✔ readr   2.1.3      ✔ forcats 0.5.2 
### ── Conflicts ─────────────────────────────────────────── tidyverse_conflicts() ──
### ✖ dplyr::filter() masks stats::filter()
### ✖ dplyr::lag()    masks stats::lag()

# You don't need to worry about the Conflicts

#################################################################
#               Diversity Tackled With R                        #
#################################################################
# Now import the data in metagenome.biom into R using the 
# import_biom() function from phyloseq
biom_metagenome <- import_biom("metagenome.biom")
# This command produces no output in the console but created a
# special class of R object which is defined by the `phloseq` 
# package and called it `biom_metagenome`. Click on the object 
# name, biom_metagenome, in the Environment pane (top right).  
# This will open a view of the object in the same pane as your script.

# A phyloseq object is a special object in R. It has five parts, called
# 'slots' which you can see listed in the object view. These are `otu_table`,
# `tax_table`, `sam_data`, `phy_tree` and `refseq`. In our case, 
# `sam_data`, `phy_tree` and `refseq` are empty. The useful data are
# in `otu_table` and `tax_table`.

# print the biom_metagenome object
biom_metagenome
# this will give the following output:
### phyloseq-class experiment-level object
### otu_table()   OTU Table:         [ 7637 taxa and 2 samples ]
### tax_table()   Taxonomy Table:    [ 7637 taxa by 7 taxonomic ranks ]

# The line starting otu_table tells us we have two samples -
# these are ERR4998593 and ERR4998600 - with a total of 7637 taxa. 
# The tax_table again tells us how many taxa we have. The seven ranks
# indicates that we have some identifications down to species level. 
# The taxonomic ranks are from the classification system of taxa from
# the most general (kingdom) to the most specific (species): kingdom/domain,
# phylum, class, order, family, genus, species.

# We can view the otu_table with:
View(biom_metagenome@otu_table)
# This will open a view of the table in the same pane as your script.

# This table has the OTU identity in the row names and the samples in the
# columns. The values in the columns are the abundance of that OTU in that
# sample.

# We can view the tax_table with:
View(biom_metagenome@tax_table)
# This gives the rank for each taxa

# k__ kingdom/domain (note a letter, two underscores)
# p__ phylum
# c__ class
# o__ order
# f__ family
# g__ genus
# s__ species

# note we have some taxa that we can't get to species level, 
# some we can't get to genus level and some we can't get to 
# family level.

# To make downstream analysis easier for us we will remove the prefix
# (e.g. f__) on each item. This contains information about the rank 
# of the assigned taxonomy, we don’t want to lose this information 
# so will and instead rename the header of each column of the DataFrame
# to contain this information.
# To remove unnecessary characters we are going to use command substring().
# This command is useful to extract or replace characters in a vector. 
# To use the command, we have to indicate the vector (x) followed by the
# first element to replace or extract (first) and the last element to 
# be replaced (last). For instance: substring (x, first, last). If a 
# last position is not used it will be set to the end of the string.
# The prefix for each item in biom_metagenome is made up of a letter 
# an two underscores, for example: o__Bacillales. In this case “Bacillales”
# starts at position 4 with an B. So to remove the unnecessary characters 
# we will use the following code:

biom_metagenome@tax_table <- substring(biom_metagenome@tax_table, 4)

# check it worked
View(biom_metagenome@tax_table)

# Let's change the names too
colnames(biom_metagenome@tax_table) <- c("Kingdom",
                                         "Phylum", 
                                         "Class", 
                                         "Order", 
                                         "Family", 
                                         "Genus", 
                                         "Species")

# check it worked
View(biom_metagenome@tax_table)


# How many OTUs are in each kingdom? We can find out by combining
# some commands. 
# We need to 
# - turn the tax_table into a data frame (a useful data structure in R)
# - group by the Kingdom column
# - summarise by counting the number of rows for each Kingdom
# This can be achieved with the following command:
biom_metagenome@tax_table %>% 
  data.frame() %>% 
  group_by(Kingdom) %>% 
  summarise(n = length(Kingdom)) 

# # A tibble: 4 × 2
# Kingdom       n
# <chr>     <int>
#   1 Archaea     339
# 2 Bacteria   7231
# 3 Eukaryota     1
# 4 Viruses      66

# Most things are bacteria.

# We can explore how many phlya we have and how many OTU there are in each
# phlya in a similar way. This time we will use View() to see the whole 
# table because it won't all print to the console. We need to 
# - turn the tax_table into a data frame (a useful data structure in R)
# - group by the Phylum column
# - summarise by counting the number of rows for each phylum
# - viewing the result
# This can be achieved with the following command:
biom_metagenome@tax_table %>% 
  data.frame() %>% 
  group_by(Phylum) %>% 
  summarise(n = length(Phylum)) %>% 
  View()
# This shows us a table with a phylum, and the number times it 
# appeared, in each row. The number of phyla is given by the 
# number of rows in this table.  By default, the table is sorted
# alphabetically by phylum. We can sort by frequency by clicking
# on the 'n' column. There are 3471 Proteobacteria and 1571
# Actinobacteria for example.

# Exercise 2
# Adapt the code to explore how many Orders we have and how many 
# OTU there are in each order.
# a) How many orders are there?
# b) What is the most common order?
# c) How many OTUs did not get down to order level?


#################################

# Plot alpha diversity
# We want to explore the bacterial diversity of our sample so 
# we will subset the bacteria

bac_biom_metagenome <- subset_taxa(biom_metagenome, Kingdom == "Bacteria")

# Now let’s look at some statistics of our bacterial metagenomes:
bac_biom_metagenome
### phyloseq-class experiment-level object
### otu_table()   OTU Table:         [ 7231 taxa and 2 samples ]
### tax_table()   Taxonomy Table:    [ 7231 taxa by 7 taxonomic ranks ]

# `phyloseq` includes a function to find the sample name and one to 
# count the number of reads in each sample.
# Find the sample names with `sample_names()`:  
sample_names(bac_biom_metagenome)
# Count the number of reads with `sample_sums()`:  
sample_sums(bac_biom_metagenome)


### "ERR4998593" "ERR4998600"
### sample_sums(bac_biom_metagenome)
### ERR4998593 ERR4998600 
###     442490     305135 

# The `summary()` function can give us an indication of species evenness  
summary(bac_biom_metagenome@otu_table)

# The median in sample ERR4998593 is 6, meaning many OTU occur six times each 
# and the maximum is very high so at least one OTU is very abundant.

# The `plot_richness()` command will give us a visual representation of 
# the diversity inside the samples (i.e. α diversity): 
plot_richness(physeq = bac_biom_metagenome,
              measures = c("Observed","Chao1","Shannon"))


# Each of these metrics can give insight of the distribution of the OTUs inside
# our samples. For example Chao1 diversity index gives more weight to singletons
# and doubletons observed in our samples, while the Shannon is a measure of 
# species richness and species evenness with more weight on richness.
# 

# Use the following to open the manual page for plot_richness
?plot_richness

# the manual page tells us how to do things like add a title, change colours,
# sort the plots etc
plot_richness(physeq = bac_biom_metagenome,
              title = "Alpha diversity indexes for both metagenomic samples",
              measures = c("Observed","Chao1","Shannon"),
              sortby = "Shannon")

#################################

# Plot beta diversity

# The β diversity between ERR4998593 and ERR498600 can be calculated with 
# the `distance()` function.

# For example Bray–Curtis dissimilarity
distance(bac_biom_metagenome, method="bray") 

# We can view our options for calculating distance with
distanceMethodList$vegdist

# e.g. looking at Jaccard distance
distance(bac_biom_metagenome, method="jaccard")

# The output of this function is a distance matrix. When we have just two 
# samples there is only one distance to calculate. If we had many 
# samples, the output would have the pairwise distances
# between all of them

#################################################################
#              Taxonomic Analysis with R                        #
#################################################################

## Reminder 
# In the last lesson, we created our phyloseq object, which contains 
# the information of our samples: `ERR4998593` and `ERR4998600`. 
# Let´s take a look again at the number of reads in our data.  
# For the whole metagenome:
biom_metagenome
sample_sums(biom_metagenome)

# Exercise 1: Repeat this for the bacterial metagenome.:

# We saw how to find out how many phyla we have and how many 
# OTU there are in each phyla by combining commands
# We 
# - turn the tax_table into a data frame (a useful data structure in R)
# - group by the Phylum column
# - summarise by counting the number of rows for each phylum
# - viewing the result
# This can be achieved with the following command:
bac_biom_metagenome@tax_table %>% 
  data.frame() %>% 
  group_by(Phylum) %>% 
  summarise(n = length(Phylum)) %>% 
  View()

#################################

# Summarise metagenomes

# `phyloseq` has a useful function that turns a phyloseq object into a dataframe.
# Since the dataframe is a standard data format in R, this makes it easier for R 
# users to apply methods they are familiar with.
# Use `psmelt()` to make a dataframe for the bacterial metagenomes
bac_meta_df <- psmelt(bac_biom_metagenome)
# Clicking on `bac_meta_df` on the Environment window will open a 
# spreadsheet-like view of it

# Now we can more easily summarise our metagenomes by sample using standard syntax.
# The following filters out all the rows with zero abundance then counts the 
# number of taxa in each phylum for each sample:

number_of_taxa <- bac_meta_df %>% 
  filter(Abundance > 0) %>% 
  group_by(Sample, Phylum) %>% 
  summarise(n = length(Abundance))
# Clicking on `number_of_taxa` on the Environment window will open a 
# spreadsheet-like view of it

# One way to visualise the phyla is with a Venn diagram. 
# The package `ggvenn` will draw one for us. It needs a data structure 
# called a list which will contain an item for each sample of the phyla
# in that sample. 
# We can see the phyla in the ERR4998600 sample with:
unique(number_of_taxa$Phylum[number_of_taxa$Sample == "ERR4998593"])


# Exercise 2: For the ERR4998600 sample

# To place the two sets of phlya in a list, we use
venn_data <- list(ERR4998593 = unique(number_of_taxa$Phylum[number_of_taxa$Sample == "ERR4998593"]),
                  ERR4998600 = unique(number_of_taxa$Phylum[number_of_taxa$Sample == "ERR4998600"]))

# And to draw the venn diagram
ggvenn(venn_data)

# The Venn diagram shows that all of the phyla are found in both samples.
# There are no phyla exclusive to either ERR4998593 or ERR4998600.
# Imagine that there were some phylas different between the two.
# Perhaps you would like to know which phyla are in ERR4998593 only? 
# The following command would print that for us:
venn_data$ERR4998593[!venn_data$ERR4998593 %in% venn_data$ERR4998600]

# Of course, in this case the result comes back with `NULL` since 
# the statement doesn't apply to any phyla.

#################################

# Visualising our data

# We summarised our metagenomes for the number of phyla in each sample 
# in `number_of_taxa`
# We can use a similar approach to examine the abundance of each of these taxa
abundance_of_taxa <- bac_meta_df %>% 
  filter(Abundance > 0) %>% 
  group_by(Sample, Phylum) %>% 
  summarise(Abundance = sum(Abundance))

# It would be nice to see the abundance by phyla as a figure to
# see if that differs between samples, even if the phyla don't.
# We can use the `ggplot()` function to visualise the breakdown 
# by Phylum in each of our two bacterial metagenomes:

abundance_of_taxa %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_col(position = "stack")


# We can see that the most abundant phyla in both ERR4998593 and ERR4998600 are 
# Proteobacteria and Actinobacteria. 

# Transformation of data
# Since our metagenomes have different sizes, it
# is imperative to convert the number of assigned read into 
# percentages (i.e. relative abundances) so as to compare them.

# We can achieve this with:

abundance_of_taxa <- abundance_of_taxa %>% 
  group_by(Sample) %>% 
  mutate(relative = Abundance/sum(Abundance) * 100)

# then plot with
abundance_of_taxa %>% 
  ggplot(aes(x = Sample, y = relative, fill = Phylum)) +
  geom_col(position = "stack")

# This bar chart still isn't ideal. It has too many
# colours to accurately identify each phyla.

# Let's try plotting the phyla as points on a scatter graph,
# with the two samples as axes.

# First we select all the columns in `abundance_of_taxa`except `Abundance`, leaving
# `Sample`, `Phylum` and `relative`. Then we pivot the data frame so that each `Sample`
# has its own column, containing the values that used to be in the `relative` column.
# One row represents one phylum.

abundance_of_taxa_wide <- abundance_of_taxa %>%
  select(-Abundance) %>%
  pivot_wider(names_from = Sample, values_from = relative)

# Then plot a scatter graph with:

abundance_of_taxa_wide %>% 
  ggplot(aes(x = ERR4998593, y = ERR4998600)) +
  geom_point() +
  ggtitle("Relative abundances of phyla")

# It looks like our data is all clustered at the low end of the scale, with only a couple
# of points further along. Let's try using a log scale to see the spread a bit better.
# We'd better add labels to our axes too, to remind us that a log scale was used

abundance_of_taxa_wide %>% 
  ggplot(aes(x = ERR4998593, y = ERR4998600)) +
  geom_point() +
  scale_x_log10() +
  xlab("Log10(ERR4998593)") +
  scale_y_log10() +
  ylab("Log10(ERR4998600)") +
  ggtitle("Relative abundances of phyla")

# That's better! We should also add a reference line too so we can see where the points
# would be if they were found in equal proportions in both samples. We can do this with 
# geom_abline() which draws a line with formula y=x

abundance_of_taxa_wide %>% 
  ggplot(aes(x = ERR4998593, y = ERR4998600)) +
  geom_point() +
  scale_x_log10() +
  xlab("Log10(ERR4998593)") +
  scale_y_log10() +
  ylab("Log10(ERR4998600)") +
  geom_abline() +
  ggtitle("Relative abundances of phyla")

# Most of the points fall closer to the ERR4998600 axis,
# telling us they were more abundant in this sample than in ERR4998593.

# Let's finish by marking the phyla that are more abundant in ERR4998593 in a different colour and 
# labelling them.
# We start by calculating the difference between the relative abundancies for each sample - this
# will tell us if each phylum is more abundant in ERR4998600 (a positive number) or in
# ERR4998593 (a negative number). Then we make a list of all the phyla which have a negative diff value

abundance_of_taxa_wide <- abundance_of_taxa_wide %>%
  mutate(diff = ERR4998600 - ERR4998593) %>%
  arrange(diff)
ERR4998593_phyla <- filter(abundance_of_taxa_wide, diff < 0)$Phylum

# Now we have a list of phyla we can plot new points on the graph that are only
# present in the list we just made, this time in a different colour.
# We can also add a label with geom_text() to tell us what phyla these points represent

abundance_of_taxa_wide %>% 
  ggplot(aes(x = ERR4998593, y = ERR4998600)) +
  geom_point() +
  geom_point(
    data = filter(abundance_of_taxa_wide, Phylum %in% ERR4998593_phyla),
    aes(), 
    col = "red") +
  geom_text(
    data = filter(abundance_of_taxa_wide, Phylum %in% ERR4998593_phyla),
    aes(label = Phylum), 
    nudge_y = -0.125,
    col = "red") +
  scale_x_log10() +
  xlab("Log10(ERR4998593)") +
  scale_y_log10() +
  ylab("Log10(ERR4998600)") +
  geom_abline() +
  ggtitle("Relative abundances of phyla")

# The two phyla which are more abundant in `ERR4998593` than in `ERR4998600` are 
# Proteobacteria and Chlamydiae. We've seen Proteobacteria before, in our bar chart, 
# where it was the most abundant phylum in both samples. We could probably have 
# guessed from that graph that Proteobacteria had a higher relative abundance in 
# one sample than the other, because the bars were very large. But we probably 
# wouldn't have been able to tell which sample Chlamydiae was more abundant in - 
# its relative proportion in both is too small to see in the bar chart. 
