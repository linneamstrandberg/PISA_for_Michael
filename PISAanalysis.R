#.rs.restartR(), some parts are unintentionally masked if codes with other libraries than
#tidyverse has been read before and therefore R session has to be restarted

library(tidyverse)
library(MSstats)

#This script is to normalize protein abundances from MS data. I check for sample outliers by an
#ecdf plot on peptide level. There is one raw file for the protein abundances and one for the peptides.
#First I convert the files to the correct format and add sample annotation (from a separate "plate order file")
#Then I remove proteins that were not detected in enough replicates. Samples are normalized based on median,
#and then q-values are calculated. Since some peptides can belong to several proteins but the normalization
#cannot handle duplicates I keep only one out of the ambiguous proteins, but add the rest after the q-values 
#have been calculated. There are probably a thousand ways to make this code shorter and more efficient,
#but my main wish is just that it works, doesn't mix up samples or something else stupid 

#I name the MS files with date, organism and "search mode" (what "search file" I used when processing the spectra)
#Then rep and tot_rep represents how many replicates of each condition and how many samples (all conditions) in total,
#since I want to make sure later on in the code that the protein was found in at least
#rep-1 samples. 

run <- "20220617"
organism <- "arabidopsis"
searchmode <- "chloroplast"
rep <- 4
tot_rep <- 8
remove_samp <- c()

#platelayout: if it reads all columns as one, change to csv2
platelayout <- read_csv(file = paste("plate/", organism, "/plate_layout_", organism, "_", run, ".csv", sep = ""))
prot <- read_tsv(paste("rawdata/", organism, "/quant_report_", organism, "_",
                      searchmode, "_", run, ".elib.proteins.txt", sep = ""))
pep <- read_tsv(paste("rawdata/", organism, "/quant_report_", organism, "_",
                      searchmode, "_", run, ".elib.peptides.txt", sep = ""))
protein_properties <- read_csv2(paste("uniprot/", organism, "/", organism, ".csv", sep = "")) %>%
  select('Entry', 'Gene names  (primary )', 'Protein names', 'Gene names  (ordered locus )')

#Check for outliers by comparing peptide intensity distribution between replicates
pep <- pep %>%
  gather(key = "Sample", value = "Intensity", 4:(ncol(prot))) %>%
  mutate(Sample = (str_sub(Sample, 50, nchar(Sample)-5))) %>%
  left_join(platelayout, by = c("Sample" = "Well")) %>%
  separate(col = Name, into = c("Treatm", "Conc", "TechReplicate"), sep = "_", remove = F) %>%
  unite("Condition", Treatm:Conc, sep = "_", remove = F)

g <- ggplot(pep, aes(log(Intensity), colour = Name)) +
      stat_ecdf(geom = "step") +
      xlab("log2(Peptide intensity)") +
      theme_bw()

rm(pep)

ggsave(filename = paste("plots/", organism, "/ecdf/ecdf_", organism, "_", searchmode, "_",
                         run, ".png", sep = ""),
        plot = g, units = "mm", width = 200, height = 80, dpi = 300)

#Annotate proteins that have ambiguous entries
prot$Ambig <- ifelse(grepl(";", prot$Protein), "Yes", "No")

prot <- prot %>% 
  separate_rows(Protein, sep = ";") %>%  
  separate(col = Protein,into = c("db", "Entry", "Entry_name"), sep = "\\|", remove = F) 

#Convert to long format. Annotate data with treatment, concentration and technical replicate
prot <- prot %>%
  gather(key = "Sample", value = "Intensity", 7:(ncol(prot)-1)) %>%
  mutate(Sample = (str_sub(Sample, 50, nchar(Sample)-5))) %>%
  left_join(platelayout, by = c("Sample" = "Well")) %>%
  separate(col = Name, into = c("Treatm", "Conc", "TechReplicate"), sep = "_", remove = F) %>%
  unite("Condition", Treatm:Conc, sep = "_", remove = F)

#annotate proteins
annot <- left_join(select(prot, Protein, Entry, Ambig, PeptideSequences),
                   protein_properties)

rm(protein_properties)

write_csv(annot, paste("annotations/", organism, "/annot_", searchmode, "_", run, ".csv", sep = ""))

amb_annot <- filter(annot, Ambig == "Yes")
write_csv(amb_annot, paste("annotations/", organism, "/amb_annot_", searchmode, "_", run, ".csv", sep = ""))

#Keep one protein if ambiguous (=peptides can be annotated to several different proteins). 
#ADD REMAINING POSSIBILITIES IN LATER STEP IN JAN'S CODE
prot <- prot %>%
  group_by(Name, Treatm) %>%
  filter(!duplicated(PeptideSequences)|n()==1)

#Filter proteins that have 0 in intensity and missing in more than one sample. 
#DO WE WANT TO FILTER PROTEINS WITH ONLY ONE PEPTIDE? - Not now, keep an eye on these
prot <- prot %>% dplyr::group_by(Condition, Protein) %>%
  dplyr::mutate(Var=sum(Intensity>0)) %>%
  filter(Var>=(rep-1)) %>% 
  select(-Var) 

prot <- prot %>% 
  group_by(Protein, Treatm) %>%
  filter(n()==tot_rep)

#Filter samples that should be removed. DOES NOT WORK NOW, FIX LATER?
# prot <- prot %>% filter(!Run %in% remove_samp)

#Prepare file for MSstats normalization with the correct column names and file format
prot <- prot %>% dplyr::rename(ProteinName = Entry, PeptideSequence = PeptideSequences, Run = Name) %>%
  mutate(PrecursorCharge = 0, FragmentIon = 0, ProductCharge = 0, IsotopeLabelType = "L", BioReplicate = 1) %>%
  select(ProteinName, PeptideSequence, PrecursorCharge,
         FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run,
         Intensity, Treatm)

prot <- prot %>% group_by(Treatm) %>% group_split()

#Normalize by medians
quant_list <- lapply(prot, function(tbl) {
  Quant <- dataProcess(
    tbl,
    logTrans = 2,
    normalization = "equalizeMedians",
    censoredInt = "0",
    cutoffCensored = "minFeature")
})

# Save Processed Data       
saveRDS(quant_list, file = paste("processeddata/", organism, "/procdata_", 
                                 organism, "_", searchmode, "_", run, ".RData", sep = ""))

#Jan's part of the code. Calculate q-values
for (n in 1:length(quant_list)) {
  
  comp <- "grouped"
  
  experiment <- substr(levels(quant_list[[n]]$ProcessedData$GROUP_ORIGINAL[1])[1], 1, nchar(levels(quant_list[[n]]$ProcessedData$GROUP_ORIGINAL[1])[1]) -2)
  
  # Create a comparison matrix for MSstats groupComparison()
  #print(levels(int_proc$ProcessedData$GROUP_ORIGINAL))
  nbr_comp <- length(levels(quant_list[[n]]$ProcessedData$GROUP_ORIGINAL)) - 1
  
  comparison_mtx <- cbind(matrix(rep(-1, nbr_comp)),
                          diag(nbr_comp))
  
  rownames(comparison_mtx) <- levels(quant_list[[n]]$ProcessedData$GROUP_ORIGINAL)[2:(nbr_comp+1)]
  
  groupComp_result <- groupComparison(contrast.matrix = comparison_mtx,
                                      data = quant_list[[n]])
  
  comparisonResult <- groupComp_result$ComparisonResult
  modelQC <- groupComp_result$ModelQC
  
  # Remove comparisons which have zero observations in any condition
  removed_comp <- comparisonResult %>% filter(!is.na(issue))
  comparisonResult <- filter(comparisonResult, is.na(issue))
  
  #Add all proteins (including the ambiguous alternatives that were removed earlier)
  stat_annot <- unique(annot) %>%
    select(Entry, PeptideSequences, `Gene names  (primary )`, Ambig) %>%
    right_join(comparisonResult, by = c("Entry" = "Protein"))
  
  stat_annot_amb <- unique(amb_annot) %>%
    select(Entry, PeptideSequences, `Gene names  (primary )`, Ambig) %>%
    left_join(stat_annot, by = "PeptideSequences") %>%
    dplyr::rename(Entry = Entry.x, `Gene names  (primary )` = `Gene names  (primary ).x`,
                  Ambig = Ambig.x) %>%
    select(-Entry.y, -`Gene names  (primary ).y`, -Ambig.y)
    
  stat <- full_join(stat_annot, stat_annot_amb) %>%
    select(-PeptideSequences)
  
  #Save
  write_csv(stat,
            paste("compresults/", organism, "/comparisonresult_", searchmode, "_", experiment, "_", comp, "_", run,
                  ".csv", sep = ""))
  write_csv(modelQC,
            paste("modelqc/", organism, "/modelqc_", searchmode, "_", experiment, "_", comp, "_", run,
                  ".csv", sep = ""))
}

#Add volcano plots
#Add sessionInfo()
