library(DSS)
library(bsseq)
library(dplyr)
library(data.table)
library(GenomeInfoDb)


# Get identified DMRs
setwd()
dss = readRDS()

dmr_coords = dss$DMR

# Get bed graphs to generate BSseq object

#MethylCount function selects the correct files corresponding to the samples within the discovery/training datasets and
#re-formats the bedGraphs for their use to build a BSseq object. This is the starting point for both DSS and DMRcate.
#The final object is a list containing: data table with samples and relative assigned bedGraph file path; bedGraph list that can be used to build a BSseq object
MethylCount = function(table, bedGraphs, datasets) {
  all_clean_sample_ids = rownames(table)
  #initialize empty vector to contain selected bedgraphs
  selected.bedGraphs = data.table()

  #singles out samples used in specific dataset and assigns the relative bedGraphs
  for (d in datasets) {
    current_plate = ifelse(d == "discovery", "Plate1",
                           ifelse(d == "p53Atypia", "Plate2",
                                  ifelse(d == "RetroValid", "Plate3",
                                         ifelse(d == "TechVal/ALQ_TWIST", "Plate4", NA))))

    # Only get the indicies of rownames which contain that plate number
    table_indicies = grep(current_plate, all_clean_sample_ids)
    sampleNum = table[table_indicies,]$sample
    bedGraph.baseName = sub("_.*", "", basename(bedGraphs[grep(d, bedGraphs)]))
    bedgraph_indicies_for_this_directory = grep(d, bedGraphs)

    bedGraph.baseName = paste0(current_plate, ".", bedGraph.baseName)
    subBedgraphs_inidicies = bedgraph_indicies_for_this_directory[match(all_clean_sample_ids, bedGraph.baseName)]
    subBedgraphs_inidicies = subBedgraphs_inidicies[!is.na(subBedgraphs_inidicies)]

    subBedgraphs = bedGraphs[subBedgraphs_inidicies]

    #add column to dataset table with the path to the relative bedGraph file
    # paths= table[grepl(d, table$plateID),][match(sampleNum,sub("_.*", "", basename(subBedgraphs))),] %>%
    #   dplyr::mutate(path2bedgraph = subBedgraphs)
    paths = table[table_indicies,] %>%
      dplyr::mutate(path2bedgraph = subBedgraphs)

    selected.bedGraphs = rbind(selected.bedGraphs, paths)
  }

  #if clause checks if the number of files is correct & if the right file has been assigned to the right sample
  if (nrow(table) != nrow(selected.bedGraphs) | all.equal(selected.bedGraphs$sample, sub("_.*", "", basename(selected.bedGraphs$path2bedgraph))) == F) {
    stop("Number of bedGraphs DOES NOT match number of samples in given dataset!", call. = FALSE)
    stop("Sample and file do not match! Double-check your input data!", call. = FALSE)

  }else {

    #format the bedGraphs in order to create a BSseq obj
    bedGraphList = lapply(lapply(selected.bedGraphs$path2bedgraph, fread), function(x) {
      dt = x %>%
        mutate(N = V5 + V6) %>% #N = total number of reads per loci
        dplyr::rename(X = V5, chr = V1, pos = V2) %>% #X = number of methylated reads per loci
        dplyr::select(chr, pos, N, X)
    })

    names(bedGraphList) = paste0(selected.bedGraphs$sample, ".", selected.bedGraphs$Case, "_", selected.bedGraphs$Pathway)

    print(head(bedGraphList))

    bedgraphObj = list(
      dataset = selected.bedGraphs,
      bedGraphList = bedGraphList
    )

    invisible(bedgraphObj)
  }
}

disc_dataset = readRDS()
datasets = c("discovery", "RetroValid", "TechVal/ALQ_TWIST", "p53Atypia")
bedGraphs = list.files(paste0("", datasets), pattern = "bedGraph", full.names = T)

bedGraph_list = MethylCount(table = disc_dataset, bedGraphs = bedGraphs, datasets = datasets)

# Now we make the BSseq object and get the smoothed values - we need to try different parameters to determine the best
# set
bs_obj = makeBSseqData(bedGraph_list$bedGraphList,
                       names(bedGraph_list$bedGraphList)) #[,]
# We need to order the object before smoothing it
bs_obj = orderBSseq(bs_obj)
seqlevels = extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="auto")
bs_obj = orderBSseq(bs_obj, seqOrder = seqlevels)

make_bsmooth_objects = function(ns = 70, maxGap = 10^8, h = 1000) {
  smoothed = BSmooth(bs_obj, ns = ns, maxGap = maxGap, h = h,
                            BPPARAM = BiocParallel::MulticoreParam(timeout = 1200, progressbar = TRUE))
  saveRDS(smoothed, paste0("BSobj_smoothed_", ns, "_", maxGap, "_", h, ".Rds"))
}

setwd()
make_bsmooth_objects()
make_bsmooth_objects(ns = 50)
make_bsmooth_objects(h = 500)