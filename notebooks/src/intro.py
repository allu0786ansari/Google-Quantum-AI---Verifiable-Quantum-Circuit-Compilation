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
# For each task the notebook tries, in order, and keeps the first that passes the
# official verifier:
#
# 1. **the public frontier** (section 1b): the strongest public circuits, built by
#    many Kaggle participants under Apache-2.0 and fully credited;
# 2. **our own circuits** (sections 2-4): invsqrt with a piecewise quadratic and
#    Booth multipliers (7.5x), dqi with tower-field arithmetic (1.53x), and a
#    chemistry pipeline;
# 3. **Google's baseline**, so the zip is always valid.
#
# Attach the competition data and ahoirg's *QCC Reproducible Circuit Assets*
# dataset. Settings are in the next cell.
