#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on Sept. 15, 2023 for running ArchR
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
## This source code is licensed under the MIT license
#######################################################################
#######################################################################

library(ArchR)
library(argparse)

parser <- ArgumentParser(description= 'The scripts to run ArchR')


parser$add_argument('--input', '-i', help= 'Fragments file folder.')
parser$add_argument('--output', '-o', help= 'Output directory.')
parser$add_argument('--genome', '-g', help= 'UCSC genome name.')
parser$add_argument('--postfix', '-f',
    help= 'Postfix of the fragments file. Defualt is .fragment.tsv.gz',
    default='.fragment.tsv.gz')
parser$add_argument('--minTSS', '-c',
    help= 'The minimum numeric transcription start site (TSS) enrichment score required for a cell to pass filtering for use in downstream analyses.',
    type= 'integer', default=4)
parser$add_argument('--minFrags', '-f',
    help= 'The minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses.',
    type= 'integer', default=1000)
parser$add_argument('--doublets', '-k',
    help= 'The number of cells neighboring a simulated doublet to be considered as putative doublets.',
    type= 'integer', default=10)
parser$add_argument('--cores', '-p', type='integer', help='Number of threads.', default=1)
parser$add_argument('--seed', '-s', type='integer', help='Random generator seed', default=1)

args<- parser$parse_args()
set.seed(args$seed)
addArchRThreads(threads=args$cores)
pf <- args$input
postfix <- args$postfix
inputFiles <- dir(pf, paste0(".*", postfix, "$"))
output <- args$output
dir.create(output, recursive=TRUE, showWarnings = FALSE)
addArchRGenome(args$genome)
# getArchRChrPrefix()
x <- read.delim((file.path(pf, inputFiles[1]), nrows=1000, header=FALSE)
xn <- names(table(x[, 1]))
if(!all(grepl("chr", xn))){
    addArchRChrPrefix(chrPrefix = FALSE)
}
ArrowFiles <- createArrowFiles(
	inputFiles = file.path(pf, inputFiles),
	sampleNames = sub(postfix, "", inputFiles),
	minTSS = args$minTSS,
	minFrags = args$minFrags,
	addTileMat = TRUE,
	addGeneScoreMat = TRUE
)
doubScores <- addDoubletScores(
	input = ArrowFiles,
	k = args$doublets, #Refers to how many cells near a "pseudo-doublet" to count.
	knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
	LSIMethod = 1
)
proj <- ArchRProject(
	ArrowFiles = ArrowFiles,
	outputDirectory = output,
	copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj <- filterDoublets(ArchRProj = proj)

p1 <- plotGroups(
	ArchRProj = proj,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "TSSEnrichment",
	plotAs = "ridges"
)
p2 <- plotGroups(
	ArchRProj = proj,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "TSSEnrichment",
	plotAs = "violin",
	alpha = 0.4,
	addBoxPlot = TRUE
)
p3 <- plotGroups(
	ArchRProj = proj,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "log10(nFrags)",
	plotAs = "ridges"
)
p4 <- plotGroups(
	ArchRProj = proj,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "log10(nFrags)",
	plotAs = "violin",
	alpha = 0.4,
	addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4,
				 name = "QC-Sample-Statistics.pdf",
				 ArchRProj = proj,
				addDOC = FALSE,
				width = 4,
				height = 4)
p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)
plotPDF(p1,p2,
				name = "QC-Sample-FragSizes-TSSProfile.pdf",
				ArchRProj = proj,
				addDOC = FALSE,
				width = 5,
				height = 5)


proj <- addIterativeLSI(
	ArchRProj = proj,
	useMatrix = "TileMatrix",
	name = "IterativeLSI")
proj <- addHarmony(
	ArchRProj = proj,
	reducedDims = "IterativeLSI",
	name = "Harmony",
	groupBy = "Sample"
)
proj <- addClusters(
	input = proj,
	reducedDims = "IterativeLSI",
	resolution = .8
	)
cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
	mat = as.matrix(cM),
	color = paletteContinuous("whiteBlue"),
	border_color = "black"
)
plotPDF(p,
				name = "sample.vs.cluster.pdf",
				ArchRProj = proj,
				addDOC = FALSE,
				width = 5,
				height = 5)

proj <- addUMAP(
	ArchRProj = proj,
	reducedDims = "IterativeLSI")
p1 <- plotEmbedding(
	ArchRProj = proj,
	colorBy = "cellColData",
	name = "Sample",
	embedding = "UMAP")
p2 <- plotEmbedding(
	ArchRProj = proj,
	colorBy = "cellColData",
	name = "Clusters",
	embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
				ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
proj <- addTSNE(
	ArchRProj = proj,
	reducedDims = "IterativeLSI",
	name = "TSNE",
	perplexity = 30
)
p1 <- plotEmbedding(
	ArchRProj = proj,
	colorBy = "cellColData",
	name = "Sample",
	embedding = "TSNE")
p2 <- plotEmbedding(
	ArchRProj = proj,
	colorBy = "cellColData",
	name = "Clusters",
	embedding = "TSNE")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf",
				ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
proj <- addUMAP(
	ArchRProj = proj,
	reducedDims = "Harmony",
	name = "UMAPHarmony",
	nNeighbors = 30,
	minDist = 0.5,
	metric = "cosine"
)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
plotPDF(p1,p2,p3,p4, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
	ArchRProj = proj,
	useMatrix = "GeneScoreMatrix",
	groupBy = "Clusters",
	bias = c("TSSEnrichment", "log10(nFrags)"),
	testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerGenes <- unique(unlist(lapply(markerList, function(.ele) head(.ele$name, n=5))))
heatmapGS <- plotMarkerHeatmap(
	seMarker = markersGS,
	cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
	labelMarkers = markerGenes,
	transpose = FALSE
)
plotPDF(
	heatmapGS,
	name = "GeneScores-Marker-Heatmap",
	width = 6, height = 9,
	ArchRProj = proj, addDOC = FALSE)
p <- plotEmbedding(
	ArchRProj = proj,
	colorBy = "GeneScoreMatrix",
	name = markerGenes,
	embedding = "UMAP",
	quantCut = c(0.01, 0.95),
	imputeWeights = NULL
)
plotPDF(plotList = p,
				name = "Plot-UMAP-Marker-Genes.pdf",
				ArchRProj = proj,
				addDOC = FALSE, width = 5, height = 5)
proj <- addImputeWeights(proj)
p <- plotEmbedding(
	ArchRProj = proj,
	colorBy = "GeneScoreMatrix",
	name = markerGenes,
	embedding = "UMAP",
	imputeWeights = getImputeWeights(proj)
)
plotPDF(plotList = p,
				name = "Plot-UMAP-Marker-Genes-Imputation.pdf",
				ArchRProj = proj,
				addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
	ArchRProj = proj,
	groupBy = "Clusters",
	geneSymbol = markerGenes,
	upstream = 50000,
	downstream = 50000
)
plotPDF(plotList = p,
				name = "Plot-Tracks-Marker-Genes.pdf",
				ArchRProj = proj,
				addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, load = FALSE)

proj <- readRDS("P7NKP/Save-ArchR-Project.rds")

pathToMacs2 <- findMacs2()
proj <- addGroupCoverages(ArchRProj = proj, groupBy="Clusters")
proj <- addReproduciblePeakSet(
	ArchRProj = proj,
	groupBy = "Clusters",
	pathToMacs2 = pathToMacs2
)
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)
markersPeaks <- getMarkerFeatures(
	ArchRProj = proj,
	useMatrix = "PeakMatrix",
	groupBy = "Clusters",
	bias = c("TSSEnrichment", "log10(nFrags)"),
	testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
head(markerList)
lengths(markerList)

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
	ArchRProj = proj,
	peakAnnotation = "Motif",
	force = TRUE
)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(
	plotVarDev,
	name = "Variable-Motif-Deviation-Scores",
	width = 5, height = 5,
	ArchRProj = proj,
	addDOC = FALSE)

saveArchRProject(ArchRProj = proj, load = FALSE)
