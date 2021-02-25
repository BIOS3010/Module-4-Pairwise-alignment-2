# Module-4---Pairwise-alignment-2
## 4.1 Manual sequence alignment of protein sequences
In this exercise, you will be working in your groups to manually generate a pairwise sequence alignment, just like we did last week. To do this, **each person in the group** draws up an alignment matrix and fills it in with the numbers and arrows between the cells. Indicate (using color or another way of higlighting) the backtracing of the optimal alignment(s). It is probably smart to find a piece of paper to draft your individual solutions. You can then either take picture of your piece of paper, or you can use the draw tool in Padlet. Feel free to use the padlet to share results with others in the group to compare your answers and verify whether they seem identical and correct.
- The Padlet you should use is here: XXX
- Each group should upload:
  1. **one** picture/drawing of an alignment matrix 
  2. The corresponding alignment(s) (use "code" formatting)
  3. Python code to generate the same alignment(s) (see below)
  4. Write "Done" at the bottom of your column, when your group is ready

Look at the provided Example to the right to see an example of how to finish the exercise

- Use the following code as inspiration. Modify the code according to your group's exercise. Use the code to check that your group's answer is correct:
```python
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.gap_score = -5
aligner.mode = 'local'

seq1 = "WLCW"
seq2= "WWGW"

alignments = aligner.align(seq1,seq2)

for alignment in alignments:
  print("Score = %.1f:" % alignment.score)
  print(alignment)
```

