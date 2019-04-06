satf_resample_create_indices = function(data, block.ids, trial.id) {
  data$row.index = 1:nrow(data)
  indices = dlply(data, block.ids, function(d) {
    indices = dlply(d, trial.id, function(d) {
      d$row.index
    })
    attr(indices, "split_type") = NULL
    attr(indices, "split_labels") = NULL
    indices
  })
  attr(indices, "split_type") = NULL
  attr(indices, "split_labels") = NULL
  indices
}

satf_resample = function(data, indices) {
  sample_indices = llply(indices, function(indices) {
    trials = sample(1:length(indices), replace=T)
    unlist(indices[trials])
  })
  data[unlist(sample_indices),]
}
