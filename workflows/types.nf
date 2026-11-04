//
// Record types shared by the local subworkflows
//
nextflow.enable.types = true

// one sample with reads; fastqs lists the runs in order: r1_run1, r2_run1, r1_run2, ... (r1 only if single-end)
record Reads {
    id: String
    single_end: Boolean
    n_runs: Integer
    fastqs: List<Path>
}

// a micom taxonomy table for one sample
record Taxonomy {
    id: String
    taxonomy: Path
}

// the growth results of one sample on one medium
record Growth {
    id: String
    medium: String
    zip: Path
    build_manifest: Path
}

// a sample that comes from the ASV table (its counts are in the table itself)
record Sample {
    id: String
}

// one diet that applies to one sample
record Diet {
    id: String
    medium: String
    diet: Path
}
