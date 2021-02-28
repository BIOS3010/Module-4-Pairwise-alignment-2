# Module-4---Pairwise-alignment-2
## 4.1 Manual sequence alignment of protein sequences
In this exercise, you will be working in your groups to manually generate a pairwise sequence alignment, just like we did last week. To do this, **each person in the group** draws up an alignment matrix and fills it in with the numbers and arrows between the cells. Indicate (using color or another way of higlighting) the backtracing of the optimal alignment(s). It is probably smart to find a piece of paper to draft your individual solutions. You can then either take picture of your piece of paper, or you can use the draw tool in Padlet. Feel free to use the padlet to share results with others in the group to compare your answers and verify whether they seem identical and correct.
- The Padlet you should use is here: https://uio.padlet.org/jonaspaulsen/m7ti9z0rhyf08wub
- Each group should upload:
  1. **one** picture/drawing of an alignment matrix 
  2. The corresponding alignment(s) (use "code" formatting)
  3. Python code to generate the same alignment(s) (see below)
  4. Write "Done" at the bottom of your column, when your group is ready

**When the results from the other groups are ready on the Padlet:**
- Discuss in the group: Are the results as expected?

To access the correct values of the substitution matrix, you can use biopython:
```python
from Bio.Align import substitution_matrices
names = substitution_matrices.load() # This is a list of the names of all available substitution matrices
mat = substitution_matrices.load("BLOSUM62")
mat["A"]["R"] # BLOSUM62 value of Alanine vs. Arginine
```


- Use the following code as inspiration. Modify the code according to your group's exercise. Use the code to check that your group's answer is correct:
```python
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.gap_score = -5

seq1 = "WLCW"
seq2 = "WWGW"

alignments = aligner.align(seq1,seq2)

for alignment in alignments:
  print("Score = %.1f:" % alignment.score)
  print(alignment)
```

## 4.2 Running BLAST using the online tool
In this group exercise, we will try to identify human homologs of a the  sequence of the RecA protein in E. coli bacteria. 

**Discuss in the group:**
- What is a homologous sequence
- Read about RecA on wikipedia: https://en.wikipedia.org/wiki/RecA 
- Discuss in the group:  Do you expect there to be a homologous sequence in human?

**Then do the following group work:**
- Open this padlet: https://uio.padlet.org/jonaspaulsen/oo4uo5iqeg9mzgo4
- Choose one in the group to share their screen.
- Find and copy the E. coli RecA protein sequence (FASTA format) from NCBI: (https://www.ncbi.nlm.nih.gov/protein)
- Go to https://blast.ncbi.nlm.nih.gov/Blast.cgi
- Select "Protein BLAST"
- Paste the copied FASTA entry into the Query sequence field
- Select Organism "Homo sapiens (taxid:9606)"
- Click "Algorithm parameters" and choose paramters according to what it says under your group in the Padlet
- Click BLAST
- Once the results are ready, describe your group's results in the Padlet using both text and/or images.

Look at the "Example" to the right of the Padlet to get an impression of how to deliver your group's answer.

**Individual work:**

## 4.3 Using affine gap penalties
Use the following code as inspiration to explore the optinal local alignments of the protein sequences `RLINLMPWVLATEYKNY` and `QFFPLMPPAPYWILATDFENY` using:
- Non-affine gap penalty: -5, matrix: BLOSUM62
- Affine gap penalties (open: -11, extend: -1), matrix: BLOSUM62
- Affine gap penalties (open: -11, extend: -1), matrix: PAM30

```python
from Bio import Align
from Bio.Align import substitution_matrices

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = 'local'
```

```diff
! Explain any differences in the optimal local alignments
! Explain the concept of affine gap penalties
! How could affine gap penalties be applied in a manual alignment?
! Advanced: Make a python script (`proteinalign.py`) that prints alignments of two input protein sequences
```

## 4.4 Running BLAST through Biopython
It is possible to run BLAST (like we did in exercise 4.2) using Biopython. This will query the online BLAST sever, but will store and keep the results of the alignment in a Python object. Note that the BLAST search will take a few miunutes, just as for the online BLAST tool.

Here is an example of how to do this, blasting a protein sequence (with blastp) towards the `nr` database at expectation value level `0.001` and showing max `10` results. The example shows a blast search of part of the RecA sequence we worked with in exercise 4.2

```python
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

seq="RLDIRRIGAVKEGENVVGSETRVKVVKNKIAAPFKQA"
result_handle = NCBIWWW.qblast(program="blastp", database="nr", expect=0.001, hitlist_size=10, sequence=seq)
blast_records = NCBIXML.parse(result_handle)
```

Just like we saw for `Bio.SeqIO` and `Bio.AlignIO` in Module 1, we have a pair of input methods, `read` and `parse`, where read is for when you have exactly one object, and parse is an iterator for when you can have lots of objects – but instead of getting e.g. a `SeqRecord` object, we get BLAST record objects.

You can use a for loop to iterate over all the results in `blast_records`:
```python
for blast_record in blast_records:
  for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        print("****Alignment****")
        print("sequence:", alignment.title)
        print("e value:", hsp.expect)
        print(hsp.align_length)
        print(hsp.query)
        print(hsp.match)
        print(hsp.sbjct)
```

```diff
! Modify the code above to do a Blast search for the entire RecA sequence
! Explain/interpret the reults
! Modify the code above to work for nucleic acids
! Blast a random nucleic acid sequence of 12 letters, and explain the results
```
