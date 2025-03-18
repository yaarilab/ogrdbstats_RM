$HOSTNAME = ""
params.outdir = 'results'  

// Process Parameters for ogrdbstats_report:
params.ogrdbstats_report.chain = params.chain

// Process Parameters for changes_names_musa:
params.changes_names_musa.chain = params.ndm_chain

// Process Parameters for changes_names_musa_1:
params.changes_names_musa_1.chain = params.ndm_chain
if (!params.airr_file){params.airr_file = ""} 
if (!params.init_run){params.init_run = ""} 
if (!params.d_genotype){params.d_genotype = ""} 
if (!params.j_genotype){params.j_genotype = ""} 
if (!params.v_genotype){params.v_genotype = ""} 
if (!params.v_germline_file){params.v_germline_file = ""} 
if (!params.musa){params.musa = ""} 
if (!params.d_germline_file){params.d_germline_file = ""} 
if (!params.j_germline_file){params.j_germline_file = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)

Channel.fromPath(params.airr_file, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_3_outputFileTSV_g_0;g_3_outputFileTSV_g_16}
Channel.fromPath(params.init_run, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_9_outputFileTSV_g_17}
Channel.fromPath(params.d_genotype, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_10_outputFileTSV_g_4}
Channel.fromPath(params.j_genotype, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_11_outputFileTSV_g_4}
Channel.fromPath(params.v_genotype, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_12_outputFileTSV_g_4}
Channel.fromPath(params.v_germline_file, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_14_germlineFastaFile_g_0;g_14_germlineFastaFile_g_19}
Channel.fromPath(params.musa, type: 'any').map{ file -> tuple(file.baseName, file) }.into{g_18_outputFileTSV_g_17;g_18_outputFileTSV_g_16}
Channel.fromPath(params.d_germline_file, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_20_germlineFastaFile_g_19}
Channel.fromPath(params.j_germline_file, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_21_germlineFastaFile_g_19}


process creat_ref_set {

input:
 set val(name1), file(v_germline_file) from g_14_germlineFastaFile_g_19
 set val(name2), file(d_germline_file) from g_20_germlineFastaFile_g_19
 set val(name3), file(j_germline_file) from g_21_germlineFastaFile_g_19

output:
 set val("reference_set"), file("${reference_set}")  into g_19_germlineFastaFile0_g_0

script:


reference_set = "reference_set_makedb.fasta"

"""
	cat ${v_germline_file} ${d_germline_file} ${j_germline_file} > ${reference_set}

"""

}


process ogrdbstats_report {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*pdf$/) "ogrdbstats_alignment/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*csv$/) "ogrdbstats_alignment/$filename"}
input:
 set val(name),file(airrFile) from g_3_outputFileTSV_g_0
 set val(name1), file(germline_file) from g_19_germlineFastaFile0_g_0
 set val(name2), file(v_germline_file) from g_14_germlineFastaFile_g_0

output:
 file "*pdf"  into g_0_outputFilePdf00
 file "*csv"  into g_0_outputFileCSV11

script:

// general params
chain = params.ogrdbstats_report.chain
outname = airrFile.name.toString().substring(0, airrFile.name.toString().indexOf("_db-pass"))

"""

germline_file_path=\$(realpath ${germline_file})

novel=""

if grep -q "_[A-Za-z][0-9]" ${v_germline_file}; then
	awk '/^>/{f=0} \$0 ~ /_[A-Za-z][0-9]/ {f=1} f' ${v_germline_file} > novel_sequences.fasta
	novel=\$(realpath novel_sequences.fasta)
	diff \$germline_file_path \$novel | grep '^<' | sed 's/^< //' > personal_germline.fasta
	germline_file_path=\$(realpath personal_germline.fasta)
	novel="--inf_file \$novel"
fi

IFS='\t' read -a var < ${airrFile}

airrfile=${airrFile}

if [[ ! "\${var[*]}" =~ "v_call_genotyped" ]]; then
    awk -F'\t' '{col=\$5;gsub("call", "call_genotyped", col); print \$0 "\t" col}' ${airrFile} > ${outname}_genotyped.tsv
    airrfile=${outname}_genotyped.tsv
fi

airrFile_path=\$(realpath \$airrfile)


run_ogrdbstats \
	\$germline_file_path \
	"Homosapiens" \
	\$airrFile_path \
	${chain} \
	\$novel 

"""

}


process changes_names_musa {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rearrangements/$filename"}
input:
 set val(name),file(airrFile) from g_3_outputFileTSV_g_16
 set val(name1),file(musa) from g_18_outputFileTSV_g_16

output:
 set val(name),file("*.tsv")  into g_16_outputFileTSV0_g_4

script:
chain = params.changes_names_musa.chain

outname = airrFile.toString()
"""
#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(stringr)
library(readr)


data <- data.table::fread("${airrFile}", data.table = F)

# Convert to data.table
setDT(data)


# Load the mapping file
mapping <- read_csv("${musa}") %>%
  select(allele, new_tag) %>%
  unique()

mapping <- mapping[grepl("${chain}",mapping[["allele"]]),]

setDT(mapping)

# Create named vector for fast lookup
replace_map <- setNames(mapping[["new_tag"]], mapping[["allele"]])


# Ensure v_call, d_call, and j_call are characters before processing
data[, v_call := sapply(v_call, function(x) {
  x <- as.character(x)  # Convert to character
  if (is.na(x)) return(NA)  # Handle NAs
  calls <- unlist(strsplit(x, ","))  # Split into components
  calls <- fcoalesce(replace_map[calls], calls)  # Replace if found
  paste(calls, collapse = ",")  # Rejoin
})]

data[, d_call := sapply(d_call, function(x) {
  x <- as.character(x)
  if (is.na(x)) return(NA)
  calls <- unlist(strsplit(x, ","))
  calls <- fcoalesce(replace_map[calls], calls)
  paste(calls, collapse = ",")
})]

data[, j_call := sapply(j_call, function(x) {
  x <- as.character(x)
  if (is.na(x)) return(NA)
  calls <- unlist(strsplit(x, ","))
  calls <- fcoalesce(replace_map[calls], calls)
  paste(calls, collapse = ",")
})]



# Write the full output file
write.table(data, sep = "\t", file = paste0("${outname}", ".tsv"), row.names = F)

"""

}


process changes_names_musa_1 {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "rearrangements/$filename"}
input:
 set val(name),file(airrFile) from g_9_outputFileTSV_g_17
 set val(name1),file(musa) from g_18_outputFileTSV_g_17

output:
 set val(name),file("*.tsv")  into g_17_outputFileTSV0_g_4

script:
chain = params.changes_names_musa_1.chain

outname = airrFile.toString()
"""
#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(stringr)
library(readr)


data <- data.table::fread("${airrFile}", data.table = F)

# Convert to data.table
setDT(data)


# Load the mapping file
mapping <- read_csv("${musa}") %>%
  select(allele, new_tag) %>%
  unique()

mapping <- mapping[grepl("${chain}",mapping[["allele"]]),]

setDT(mapping)

# Create named vector for fast lookup
replace_map <- setNames(mapping[["new_tag"]], mapping[["allele"]])


# Ensure v_call, d_call, and j_call are characters before processing
data[, v_call := sapply(v_call, function(x) {
  x <- as.character(x)  # Convert to character
  if (is.na(x)) return(NA)  # Handle NAs
  calls <- unlist(strsplit(x, ","))  # Split into components
  calls <- fcoalesce(replace_map[calls], calls)  # Replace if found
  paste(calls, collapse = ",")  # Rejoin
})]

data[, d_call := sapply(d_call, function(x) {
  x <- as.character(x)
  if (is.na(x)) return(NA)
  calls <- unlist(strsplit(x, ","))
  calls <- fcoalesce(replace_map[calls], calls)
  paste(calls, collapse = ",")
})]

data[, j_call := sapply(j_call, function(x) {
  x <- as.character(x)
  if (is.na(x)) return(NA)
  calls <- unlist(strsplit(x, ","))
  calls <- fcoalesce(replace_map[calls], calls)
  paste(calls, collapse = ",")
})]



# Write the full output file
write.table(data, sep = "\t", file = paste0("${outname}", ".tsv"), row.names = F)

"""

}

g_10_outputFileTSV_g_4= g_10_outputFileTSV_g_4.ifEmpty([""]) 

def defaultIfInexistent(varName){
    try{
    	println binding.hasVariable(varName)
        varName.toString()
        println varName()
        return varName
    }catch(ex){
        return "check"//file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
    }
}

def bindingVar(varName) {
    def optVar = binding.hasVariable(varName)//binding.variables.get(varName)
    println optVar
    if(optVar) {
    	println "pass"
        println optVar
        //will only run for global var
    }
    println "fail"
    optVar
}
process VDJbase_genotype_report {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${outname}_Final_genotype.tsv$/) "genotype_report/$filename"}
input:
 set val(name1),file(initial_run) from g_17_outputFileTSV0_g_4
 set val(name2),file(personal_run) from g_16_outputFileTSV0_g_4
 set val(name3),file(v_genotype) from g_12_outputFileTSV_g_4
 set val(name4),file(d_genotype) from g_10_outputFileTSV_g_4
 set val(name5),file(j_genotype) from g_11_outputFileTSV_g_4

output:
 set val(outname),file("${outname}_Final_genotype.tsv")  into g_4_outputFileTSV00

script:

outname = initial_run.name.substring(0, initial_run.name.indexOf("_Second_Alignment_db-pass"))

"""
#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(alakazam)

# the function get the alleles calls frequencies
getFreq <- function(data, call = "v_call"){
	# get the single assignment frequency of the alleles
	table(grep(",", data[[call]][data[[call]]!=""], invert = T, value = T))
}

addFreqInfo <- function(tab, gene, alleles){
	paste0(tab[paste0(gene, "*", unlist(strsplit(alleles, ',')))], collapse = ";")
}

## read selected data columns

data_initial_run <- fread("${initial_run}", data.table = FALSE, select = c("sequence_id", "v_call", "d_call", "j_call"))
data_genotyped <- fread("${personal_run}", data.table = FALSE, select = c("sequence_id", "v_call", "d_call", "j_call"))

## make sure that both datasets have the same sequences. 
data_initial_run <- data_initial_run[data_initial_run[["sequence_id"]] %in% data_genotyped[["sequence_id"]],]
data_genotyped <- data_genotyped[data_genotyped[["sequence_id"]] %in% data_initial_run[["sequence_id"]],]
data_initial_run <- data_initial_run[order(data_initial_run[["sequence_id"]]), ]
data_genotyped <- data_genotyped[order(data_genotyped[["sequence_id"]]), ]

non_match_v <- which(data_initial_run[["v_call"]]!=data_genotyped[["v_call"]])

data_initial_run[["v_call"]][non_match_v] <- data_genotyped[["v_call"]][non_match_v]
    

# for the v_calls
print("v_call_fractions")
tab_freq_v <- getFreq(data_genotyped, call = "v_call")
tab_clone_v <- getFreq(data_initial_run, call = "v_call")
# keep just alleles that passed the genotype
tab_clone_v <- tab_clone_v[names(tab_freq_v)]
# read the genotype table
genoV <- fread("${v_genotype}", data.table = FALSE, colClasses = "character")
# add information to the genotype table
genoV <-
  genoV %>% dplyr::group_by(gene) %>% dplyr::mutate(
    Freq_by_Clone = addFreqInfo(tab_clone_v, gene, GENOTYPED_ALLELES),
    Freq_by_Seq = addFreqInfo(tab_freq_v, gene, GENOTYPED_ALLELES)
  )


# for the j_calls
print("j_call_fractions")
tab_freq_j <- getFreq(data_genotyped, call = "j_call")
tab_clone_j <- getFreq(data_initial_run, call = "j_call")
# keep just alleles that passed the genotype
tab_clone_j <- tab_clone_j[names(tab_freq_j)]
# read the genotype table
genoJ <- fread("${j_genotype}", data.table = FALSE, colClasses = "character")
# add information to the genotype table
genoJ <-
  genoJ %>% dplyr::group_by(gene) %>% dplyr::mutate(
    Freq_by_Clone = addFreqInfo(tab_clone_j, gene, GENOTYPED_ALLELES),
    Freq_by_Seq = addFreqInfo(tab_freq_j, gene, GENOTYPED_ALLELES)
  )
  
# for the d_calls; first check if the genotype file for d exists
# if("${d_genotype}"=="*tsv")
if (endsWith("${d_genotype}", ".tsv")){
	# for the d_calls
	print("d_call_fractions")
	tab_freq_d <- getFreq(data_genotyped, call = "d_call")
	tab_clone_d <- getFreq(data_initial_run, call = "d_call")
	# keep just alleles that passed the genotype
	tab_clone_d <- tab_clone_d[names(tab_freq_d)]
	# read the genotype table
	genoD <- fread("${d_genotype}", data.table = FALSE, colClasses = "character")
	# add information to the genotype table
	print(tab_clone_d)
	print(tab_freq_d)
	print(genoD)
	genoD <-
	  genoD %>% dplyr::group_by(gene) %>% dplyr::mutate(
	    Freq_by_Clone = addFreqInfo(tab_clone_d, gene, GENOTYPED_ALLELES),
	    Freq_by_Seq = addFreqInfo(tab_freq_d, gene, GENOTYPED_ALLELES)
	  )
	  
	genos <- plyr::rbind.fill(genoV, genoD, genoJ)
}else{
	genos <- plyr::rbind.fill(genoV, genoJ)
}

genos[["Freq_by_Clone"]] <- gsub("NA", "0", genos[["Freq_by_Clone"]])
genos[["Freq_by_Seq"]] <- gsub("NA", "0", genos[["Freq_by_Seq"]])

# rename the genotyped_allele columns
new_genotyped_allele_name = "GENOTYPED_ALLELES"
col_loc = which(names(genos)=='GENOTYPED_ALLELES')
names(genos)[col_loc] = new_genotyped_allele_name


# write the report
print("Writing Genotype Report")
write.table(genos, file = paste0("${outname}","_Final_genotype.tsv"), row.names = F, sep = "\t")
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
