library(RJSONIO)
library(googlesheets4)

# Set command line arguments
args <- commandArgs(TRUE)

# Path to Google Sheet. Link should be made public
gsheet_path <- "https://docs.google.com/spreadsheets/d/1cUccrex6xkkyhQkaN__3xmjn6BVu5BK_BsrnSBqhGd8/edit?usp=sharing"

#gsheet_path <- args[1]

# Out path to json file for loading into IGV
out_path <- "browser/neuro-dev-igv.json"
#out_path <- args[2]

# Specify reference genome.
genome <- "hg19"

#genome <- args[3]
genome <- as.character(genome)

supported_genomes <- c("hg19", "mm10")

genome_pass <- genome %in% supported_genomes

message("##############")

if (genome_pass == TRUE){
    yes = message("Reading Google Sheet...")
} else {
    message(paste("Specified genome incorrect. Must be one of the following:",
                  paste(supported_genomes, collapse = ", ")))
}

stopifnot(genome_pass)

message("##############")


# Read the google sheet data
gs4_deauth() # ensures google does not need authentication to read sheet
mdat <- read_sheet(gsheet_path, col_types = "cccclccccilciiic")
mdat <- data.frame(mdat)


# Subset data based on T/F for include in browser column
mdat <- mdat[mdat$Include_in_browser, ]

# Remove rows with no links.
mdat <- mdat[!is.na(mdat$url), ]

message("##############")
message("Creating hg19 IGV session file...")
message("##############")

# setup for human hg19
hg19_seed <- c('{"reference": {
		"id": "hg19",
		"name": "Human (CRCh37/hg19)",
		"fastaURL": "https://s3.dualstack.us-east-1.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta",
		"indexURL": "https://s3.dualstack.us-east-1.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/hg19.fasta.fai",
		"cytobandURL": "https://s3.dualstack.us-east-1.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg19/cytoBand.txt"
	},
	"tracks": [
		{
			"type": "sequence",
			"order": -1.797693e+308,
			"noSpinner": true
		},
		{
			"name": "Refseq Genes",
			"format": "refgene",
			"url": "https://s3.dualstack.us-east-1.amazonaws.com/igv.org.genomes/hg19/refGene.sorted.txt.gz",
			"indexURL": "https://s3.dualstack.us-east-1.amazonaws.com/igv.org.genomes/hg19/refGene.sorted.txt.gz.tbi",
			"visibilityWindow": -1,
			"removable": false,
			"order": 1,
			"noSpinner": true,
			"filename": "refGene.sorted.txt.gz",
			"sourceType": "file",
			"type": "annotation",
			"maxRows": 500,
			"filterTypes": [
				"chromosome",
				"gene"
			]
		}
	]
}')

mm10_seed <- c('{
	"version": "2.7.2",
	"reference": {
		"id": "mm10",
		"name": "Mouse (GRCm38/mm10)",
		"fastaURL": "https://s3.dualstack.us-east-1.amazonaws.com/igv.broadinstitute.org/genomes/seq/mm10/mm10.fa",
		"indexURL": "https://s3.dualstack.us-east-1.amazonaws.com/igv.broadinstitute.org/genomes/seq/mm10/mm10.fa.fai",
		"cytobandURL": "https://s3.dualstack.us-east-1.amazonaws.com/igv.broadinstitute.org/annotations/mm10/cytoBandIdeo.txt.gz",
		"order": 1000000
	},
	"tracks": [
		{
			"name": "Refseq Genes",
			"format": "refgene",
			"url": "https://s3.dualstack.us-east-1.amazonaws.com/igv.org.genomes/mm10/ncbiRefSeq.txt.gz",
			"indexed": false,
			"order": 1000000,
			"removable": false,
			"visibilityWindow": -1,
			"supportsWholeGenome": false
		}
	]
}')


# Set json seed based on argument 2

seed <- NULL

if (genome == "hg19"){
    seed <- hg19_seed
} else if (genome == "mm10"){
    seed <- mm10_seed
}


start_json <-RJSONIO::fromJSON(seed)
tracks <- start_json$tracks

# Function to make tracks info
make_json <- function(x, filename=out_path, sourceType="file",
                      priority_offset=0, noSpinner=TRUE){

    new_track <- mdat[x, ]

    # Set track name
    name <- as.character(new_track$track_name)

    # Set file path
    url <- as.character(new_track$url)

    # Set track order
    order <- new_track$order + priority_offset

    # Construct the output
    new_track_info <- list(name=name, filename=filename, format=new_track$format,
                           url=url, sourceType=sourceType,
                           type=new_track$type, noSpinner=noSpinner,
                           color=new_track$color,
                           height=new_track$height, autoscale=new_track$autoscale,
                           autoscaleGroup=new_track$autoscaleGroup,
                           min=new_track$min, max=new_track$max,
                           order=new_track$order+priority_offset)
    return(new_track_info)

}

### Make tracks
data_tracks <- lapply(X = 1:nrow(mdat), FUN = make_json)
all_tracks <- c(tracks, data_tracks)
start_json$tracks <- all_tracks

#save to json format
exportJSON <- RJSONIO::toJSON(start_json, pretty = TRUE)
write(exportJSON, out_path)

message("##############")
message("Done!")
message("##############")
message(paste("IGV browser session file saved to", out_path))
message("##############")
