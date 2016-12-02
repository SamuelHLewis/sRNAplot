# sRNAplot
## Purpose
A tool to visualize various small RNA characteristics
## Requirements
Written in python3.

Requires:

[WebLogo 3](http://weblogo.berkeley.edu/)

[signature.py](http://52.203.54.162/GEDlab/?page_id=730) by Christophe Antoniewski (must be installed as ~/bin/signature.py)

## Usage
Basic usage is:
```bash
python3 sRNAplot.py -s samplename -i input.bam
```
sRNAplot accepts the following additional arguments (all of which have defaults already set):

-m (minimum sRNA length, default=17)

-M (maximum sRNA length, default=35)

--seqlogo (plot separate seqlogos for sense and antisense sRNAs)

--pingpong (plot Ping-Pong overlap signature)

--unique (plot unique reads only i.e. collapse identical reads)

-h display help message
