library("phylotools")
library(data.table)
library(dplyr)
library(stringr)
library(Biostrings)

# Argument parser
parser <- ArgumentParser(description = 'Combine IG/TCR annotations from IgDetective & Digger.')

parser$add_argument("--base_dir", required=TRUE,
                    help="Base directory path")
parser$add_argument("--locus", required=TRUE,
                    help="locus")

args <- parser$parse_args()

# Assign arguments
base_dir<- args$base_dir
locus <- args$locus

igd_path<-paste0(base_dir,"combined_genes_",locus,".txt")
igd<-fread(igd_path)

digger_pathfw<-paste0(base_dir,locus,"_digger_f.tsv")
digger_pathrv<-paste0(base_dir,locus,"_digger_r.tsv")

if (file.exists(digger_pathfw) & file.exists(digger_pathrv)) {
  digger <- rbind(fread(digger_pathfw), fread(digger_pathrv))
} else if (file.exists(digger_pathfw)) {
  digger <- fread(digger_pathfw)
} else if (file.exists(digger_pathrv)) {
  digger <- fread(digger_pathrv)
} else {
  stop("Neither digger_pathfw nor digger_pathrv exists.")
}


digger <- digger %>%
  mutate(contig_start = as.numeric(str_extract(contig, "(?<=start_pos)\\d+")),
         Contig = str_extract(contig, "(?<=contig)[^s]+"),
         Productive = case_when(
           functional == "Functional" ~ TRUE,
           functional == "ORF" ~ TRUE,
           functional == "pseudo" ~ FALSE
         ))


digger$start_pos<-digger$contig_start+digger$start-1
digger$end_pos<-digger$contig_start+digger$end
digger$Locus<-locus
igd_regions <- as.data.frame(igd %>%
                               mutate(source = "IgDetective",
                                      End = Pos + nchar(Sequence)))

digger_regions <- as.data.frame(digger %>%
                                  mutate(source = "Digger",
                                         GeneType = gsub(locus,"",gene_type),
                                  ) %>%
                                  select(GeneType,Contig, Pos = start_pos, 
                                         Strand = sense, Sequence=seq,
                                         Productive, Locus,source,End = end_pos))


overlap_combinations <- expand.grid(
  igd_index = seq_len(nrow(igd_regions)),
  digger_index = seq_len(nrow(digger_regions))
)

overlap_filtered <- overlap_combinations %>%
  filter(
    igd_regions$Pos[igd_index] <= digger_regions$End[digger_index] &
      igd_regions$End[igd_index] >= digger_regions$Pos[digger_index]
  )

overlap_igd <- igd_regions[unique(overlap_filtered$igd_index), ]
overlap_digger <- digger_regions[unique(overlap_filtered$digger_index), ]

overlap_all <- bind_rows(
  overlap_igd %>% mutate(source = "Both"),
  overlap_digger %>% mutate(source = "Both")
)

igd_overlap_indices <- unique(overlap_filtered$igd_index)
digger_overlap_indices <- unique(overlap_filtered$digger_index)

# IGDetective only: rows in igd_regions not in overlap_filtered
igd_only_regions <- igd_regions %>%
  dplyr::slice(-igd_overlap_indices)

# Digger only: rows in digger_regions not in overlap_filtered
digger_only_regions <- digger_regions %>%
  dplyr::slice(-digger_overlap_indices)

both_regions<- igd_regions %>%
  dplyr::slice(igd_overlap_indices)%>% mutate(source = "IgDetective+Digger")

digger_only_regions$Contig<- toupper(sub("(.*)(\\d)$", "\\1.\\2", digger_only_regions$Contig))

all<-rbind(both_regions,digger_only_regions,igd_only_regions)

productive_all<-all%>% 
  mutate(Productive = TRUE)

write_csv(productive_all,paste0(base_dir,locus,"_combined_productive.csv"))
write_csv(all,paste0(base_dir,locus,"_combined.csv"))

