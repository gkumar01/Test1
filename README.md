# Test1

#Requirements
1. install bedtools:
sudo apt-get install bedtools

2. Python package installation:
pip install numpy
pip install pysam
pip install pybedtools
pip install pandas


#3. Input files:
1. TEST1.bam : aligned read bam file
2. TEST1_region.bed : target sites 

# Running python script
python test_bam.py -i TEST1.bam -t TEST1_region.bed