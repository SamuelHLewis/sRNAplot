# sRNAplot
## Purpose
A tool to visualize various small RNA characteristics
## Requirements
Written in python3.

Requires:

[WebLogo 3](http://weblogo.berkeley.edu/)

signature.py by [Christophe Antoniewski](http://52.203.54.162/GEDlab/) - to install, simply copy the `signature.py` file (above) into your `~/bin/` directory (**NOTE**: it **must** be installed in this directory).

## Usage
Basic usage is:
```bash
sRNAplot.py -s samplename -i input.bam
```
sRNAplot accepts the following additional arguments (all of which have defaults already set):

-m (minimum sRNA length, default=17)

-M (maximum sRNA length, default=35)

--seqlogo (plot separate seqlogos for sense and antisense sRNAs)

--pingpong (plot Ping-Pong overlap signature)

--distribution (plot sRNA length distribution)

--unique (plot unique reads only i.e. collapse identical reads)

-h (display help message)
