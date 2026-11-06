# %% [markdown]
# # Verifiable Quantum Circuit Compilation: submission notebook
#
# This notebook builds all three circuits from source, checks them, and writes
# `submission.zip` in the layout the competition expects:
#
# ```
# invsqrt.kmb
# dqi.kmb
# dfthcsos/angle_k_0.kmb ... angle_k_<K-1>.kmb, submission.npz
# ```
#
# | Task | Version in this notebook |
# |---|---|
# | invsqrt | **Improved**: piecewise quadratic replaces the Newton step (section 2) |
# | dqi | Baseline reference decoder (improvement in progress) |
# | dfthcsos | Baseline sample submission (improvement in progress) |
#
# Every task falls back to its baseline circuit if a check fails, so the zip is
# always valid. Settings are in the next cell.
