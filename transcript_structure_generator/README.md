# Setup

To install the Python virtual environment, run

```
conda env create --file environment.yml
conda activate transcript-structure-generator
```

# Usage

Input:
- Csv-formatted file ("ID,Count") with counts for individual transcripts
- Probability of intron inclusion (float in range [0,1])
- gtf-formatted file with exon coordinates of the transcripts included in the csv file

Output:
- gtf-formatted file containing generated intron/exon structures per transcript
- csv-formatted file ("NewTranscriptID,ID,Count") with
	- id of generated transcript
	- id of original transcript (without intron inclusions)
	- count