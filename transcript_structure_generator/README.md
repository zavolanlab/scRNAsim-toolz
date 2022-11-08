# Setup

To install the Python virtual environment, run

```
conda env create --file environment.yml
conda activate transcript-structure-generator
```

# Usage

Input:
- csv-formatted file ("ID,Count") with counts for individual transcripts
- probability of intron inclusion (float in range [0,1])
- gtf-formatted file with exon coordinates of the transcripts included in the csv file

Output:
- gtf-formatted file containing generated intron/exon structures per transcript
- csv-formatted file ("NewTranscriptID,ID,Count") with
	- id of generated transcript
	- id of original transcript (without intron inclusions)
	- count

To install package, run

```
pip install .
```

Afterwards, it can be imported using

```python
import tsg
```

To generate the sampled transcripts, run

```
transcript-generator --transcripts <transcripts_file> --annotation <annotations_file> --prob_inclusion=<probability_inclusion>
```

where the transcripts file should be csv-formatted, the annotation file gtf-formatted and the inclusion probability for introns a float in the range [0,1].


# Development

To perform all tests, make sure your environment corresponds to the `environment.yml` file and run

```
pytest tests
```