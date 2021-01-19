filter_with_cdhit <- function(sequences, threshold, word_length = 2, cdhit_path = "./third-party") {
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  cdhit <- paste0(cdhit_path, "/cdhit -i ", input,  " -o ", output, " -c ", threshold, " -n ", word_length)
  write_fasta(sequences, input)
  system(cdhit)
  res <- read_fasta(output)
  file.remove(input, output, paste0(output, ".clstr"))
  res
}
