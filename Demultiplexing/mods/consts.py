"""Constant numbers used in remove-background."""

# All available methods
METHODS = ["popscle", "souporcell", "DoubletDetection", "DoubletFinder", "scDblFinder", "scds", "Scrublet"]

# Demultiplexing methods.
DEMULTIPLEX = ["popscle", "souporcell"]

# Default doublet detection methods for single-cell or single-nucleus.
SINGLE_CELL = ["DoubletDetection", "scds", "Scrublet"]
SINGLE_NUCLEUS = ["DoubletFinder", "scDblFinder"]