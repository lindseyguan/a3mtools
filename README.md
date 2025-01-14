# a3mtools

Tools for working with a3m files. Designed for generating input to structure prediction tools like alphafold.

## Installation

```bash
git clone https://github.com/jacksonh1/a3mtools.git
cd a3mtools
pip install .
```

## What these tools do:

The current implementation of the tools are designed to handle MSAs retrieved using the colabfold MMseqs2 search tool. There are some strange quirks about those files specifically that I tried to account for. Mainly query sequence names are denoted by 3 digit integers starting at 101. So if an MSA only has 1 query sequence it will be named 101. If there are 2 query sequences they will be named 101 and 102 etc. This comes into play when combining MSAs for predicting protein complexes. Additionally, when combining 2 MSAs, the query sequences are combined into a single sequence. This is required to predict the complex structure. The query sequences are also added back to the MSA as individual sequences. **MSAs are combined in unpaired format**.<br>
<br>
Here's an example (with no homologous sequences):
<br>
MSA1 a3m file:
```
#9	1
>101
ABCDEFGHI
```
MSA2 a3m file:
```
#17	1
>101
JKLMNOPQRSTUVWXYZ
```
MSA1 + MSA2 a3m file:
```
#9,17   1,1
>101    102
ABCDEFGHIJKLMNOPQRSTUVWXYZ
>101
ABCDEFGHI-----------------
>102
---------JKLMNOPQRSTUVWXYZ
```
The MSA1 + MSA2 a3m file could be used as input to many structure prediction tools to predict the structure of the complex formed by the 2 query sequences. <br>

I don't know if maintaining this specific 101, 102, 103 ... naming scheme is strictly necessary, but we'll stick with it for now. <br>

## Usage
There will eventually be a more complete guide. But for now, you can install the package and run the following code to see some of the basic functionality. <br>


```python
import a3mtools
import a3mtools.examples as examples

# import an a3m file
msa = a3mtools.MSAa3m.from_a3m_file(examples.example_a3m_file1)
print(msa)

# slicing the alignment
print(msa[2:5])

# concatenating alignments
msa2 = a3mtools.MSAa3m.from_a3m_file(examples.example_a3m_file2)
print(msa2)
print(msa + msa2)
print(msa + msa2 + msa)

# saving the alignment to a file
complex_msa = msa + msa2
complex_msa.save("example_complex.a3m")
```

## future features:
---
- [ ] documentation on readthedocs
  - [ ] examples and code
- [ ] convert between fasta and a3m and back
- [ ] allow for more generic naming of query sequence
- [ ] add better test functions
- [ ] an option for combining MSAs in paired format

