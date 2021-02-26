# _De novo_ microsatellite loci identification pipeline

## Description
__Identification of microsatellite loci together with corresponding primer sequences and validating specificity and potential for multiplexing.__

Microsatellite loci are identified with MISA within a sequence set (e.g. *de novo* assembly). Afterwards the results are filtered for non-compound loci with a maximum of 20 repetitions. The Primer3 interface modules of MISA where customized to be compatible with the used release version (`libprimer3 release 2.3.7`). The script `p3_in.pl` was adjusted to `my_p3_in.pl` to adapt and extend the options for Primer3. Afterwards `primer3_core` can be run with the corresponding output. The script `p3_out.pl` was rewritten as `my_p3_out.pl` to return information on primer sequences and melting temperatures, fragment size, motif and its length from the MISA and Primer3 output files in tabular format.
To ensure accuracy the primer sequences are converted to fasta format and are aligned to the sequence set they were extracted from via `blastn`. The resulting blast hits are filtered for error free primer alignments over the complete primer. Primer pairs are kept only if both sequences survive filtering.
To ensure specificity and potential for multiplexing the remaining blast hits are converted to `bed` format and checked if there are no primer pairs binding within 3kb up- and downstream as well as there is a larger distance than 3kb to the sequence edge.

## Dependencies
- [MISA](https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_22092015.zip)
- [Primer3](https://github.com/primer3-org/primer3)
- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)

## Input
A set of sequences to screen for microsatellite loci. Be sure that the sequences potentially fulfill the requirements of your parameters, which is distance of a locus to a sequence edge.

## Usage
There is no script (yet) to execute the whole pipeline at once. You can follow the commands below and adjust parameters manually for now.

```
#run MISA on the sequence set
#set parameters can be found in misa.ini
misa.pl contigs_ge100.fasta > misa.log 2> misa.err

#filter MISA output for non-compound loci with a maximum of 20 repetitions
awk -F'[\t)]' '$3!~/^c/&&$5<=20' contigs_ge100.fasta.misa > contigs_ge100.fasta.misa_filter

#create primer3 input file
my_p3_in.pl contigs_ge100.fasta.misa_filter

#run primer3
primer3_core contigs_ge100.fasta.p3in > contigs_ge100.fasta.p3out

#extract information from primer3 output
my_p3_out.pl contigs_ge100.fasta.misa_filter contigs_ge100.fasta.p3out > contigs_ge100.fasta.p3out.txt

#convert primer sequences into fasta format
awk '$1!="ID"{print ">"$1"_F\n"$2"\n>"$1"_R\n"$4}' contigs_ge100.fasta.p3out.txt > contigs_ge100.fasta.p3out.primer.fasta

#blast primer sequences against the sequence set
makeblastdb -in ./contigs_ge100.fasta -out ./contigs_ge100 -dbtype nucl > mbdb.log 2> mbdb.err
blastn -db ./contigs_ge100 -query ./contigs_ge100.fasta.p3out.primer.fasta -out contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100 -outfmt '6 std qlen' -evalue 1000 -word_size 7 -dust no -num_threads 80 > blastn.log 2> blastn.err

#filter for 100% identity matches over the full length of the query
awk '$3==100&&$4==$13' contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100 > contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches

#extract specific primer pair alignments
for i in `cut -f1 contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches | sed 's/_[F,R]$//' | sort | uniq -c | awk '$1==2{print $2}'`; do grep $i contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches; done > contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1

#extract one primer pair per locus
#one can limit the returned pairs in primer3 as well but you never know if the first returned pairs are always specific in your sequence set
for i in `cut -f1 contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1 | sed 's/_[0-4]_[F,R]$//' | sort -u`; do grep $i contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1 | head -n 2; done > contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1.1st-pair

#convert the blast hits to bed format
awk '{if($9>$10){print $2"\t"$10-1"\t"$9"\t"$1}else{print $2"\t"$9-1"\t"$10"\t"$1}}' contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1.1st-pair > contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1.1st-pair.bed
#sort and verify bed format
bedtools sort -i contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1.1st-pair.bed > contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1.1st-pair.sort.bed

#check multiplexing potential
#no other primer pair should bind within 3kb up- and downstream
bedtools merge -i contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1.1st-pair.sort.bed -d 3000 -c 4 -o collapse | awk -F'[\t,]' 'NF==5{gsub("_F","",$4);print $4}' > contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1.1st-pair.sort.bed.candidates
#no locus should be closer to a sequence edge than 3kb
awk -F'[_\t]' '$11>3000&&$4-$12>3000' contigs_ge100.fasta.misa_filter > contigs_ge100.fasta.misa_filter.candidates

#get a list of ids that fulfill all filtering criteria
for i in `awk '{print $1"_"$2}' contigs_ge100.fasta.misa_filter.candidates`; do grep $i contigs_ge100.fasta.p3out.primer.fasta_vs_contigs_ge100_100perc-matches.spec1.1st-pair.sort.bed.candidates; done > consensus.candidates

#return the final output
grep -f consensus.candidates contigs_ge100.fasta.p3out.txt | sort -nk8 > contigs_ge100.fasta.p3out.consensus.candidates.txt
```

### Output
The final output is a `tsv` file containing:
- __Locus ID:__ Combined ID of `fasta` ID from sequence set, locus ID from MISA and primer pair ID from Primer3.
- __Forward sequence:__ Nucleotide sequence of the forward primer.
- __Forward Tm:__ Melting temperature of the forward primer.
- __Reverse sequence:__ Nucleotide sequence of the reverse primer.
- __Reverse Tm:__ Melting temperature of the reverse primer.
- __Length:__ Fragment size of the locus.
- __Motif:__ Nucleotide sequence and repetitions of the identified locus.
- __Motif length:__ Length of the microsatellite sequence in bp.

## Citation
If you use this pipeline please cite:

- Schröder O, Cavanaugh KK, Schneider JV, Schell T, Bonada N, Seifert L, Pauls SU (2021). Genetic data support local persistence in multiple glacial refugia in the montane net‐winged midge Liponeura cinerascens cinerascens (diptera, blephariceridae). Freshwater Biology. <https://doi.org/10.1111/fwb.13682>

Additional to the dependencies:

- MISA:  
Beier S, Thiel T, Münch T, Scholz U, Mascher M (2017). MISA-web: a web server for microsatellite prediction. _Bioinformatics_, 33(16):2583–2585. <https://doi.org/10.1093/bioinformatics/btx198>
- Primer3:  
Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC et al. (2012) Primer3—new capabilities and interfaces. _Nucleic Acids Research_, 40(15):e115. <https://doi.org/10.1093/nar/gks596>
- blast:  
Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J et al. (2009). BLAST+: architecture and applications. _BMC Bioinformatics_, 10(1):421. <https://dx.doi.org/10.1186%2F1471-2105-10-421>
- bedtools:   
Quinlan AR, Hall IM (2010). BEDTools: a flexible suite of utilities for comparing genomic features. _Bioinformatics_, 26(6):841–842. <https://doi.org/10.1093/bioinformatics/btq033>
