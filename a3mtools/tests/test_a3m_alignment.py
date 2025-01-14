# %%

import a3m_tools
import a3m_tools.config as config

# ==============================================================================
# // EXAMPLE USAGE
# ==============================================================================

# import an a3m file
msa = a3m_tools.MSAa3m.from_a3m_file(config.example_a3m_files[0])
print(msa)

# slicing the alignment
print(msa[2:5])

# concatenating alignments
msa2 = a3m_tools.MSAa3m.from_a3m_file(config.example_a3m_files[1])
print(msa2)
print(msa + msa2)
print(msa + msa2 + msa)

# saving the alignment
msa[2:5].save("test_alignment_1_sliced.a3m")
