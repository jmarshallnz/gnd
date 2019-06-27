## Running error correction model

This assumes we're running the error correction where there's no ground-truth. i.e. no controls with known
concentrations that we can use to estimate error rates.

### 1. Run `create_distance_matrix.R`

You may need to alter this as it's setup for the 2 farm stuff. Basically the idea is this creates the distance
matrix between pairs. The `read_fasta_path` method might be better than `read_fasta` in some cases.

It also creates a data.frame with closest pairs.

### 2. Run `read_errors.R`
