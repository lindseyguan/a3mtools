# %%

import a3m_tools

msa = a3m_tools.MSAa3m.from_a3m_file("./test_alignment_1.a3m")
msa2 = a3m_tools.MSAa3m.from_a3m_file("./test_alignment_2.a3m")

# %%
# ==============================================================================
# // EXAMPLE USAGE
# ==============================================================================

# import an a3m file
msa = a3m_tools.MSAa3m.from_a3m_file("./test_alignment_1.a3m")
print(msa)

# slicing the alignment
print(msa[2:5])

# concatenating alignments
msa2 = a3m_tools.MSAa3m.from_a3m_file("./test_alignment_2.a3m")
print(msa2)
print(msa + msa2)
print(msa + msa2 + msa)

# saving the alignment
msa[2:5].save("test_alignment_1_sliced.a3m")
