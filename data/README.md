# Data

Data is stored elsewhere (TODO: add accession numbers when data is deposited.)

File structure should be:
data
|── MhiMlo 
|        |── Mhi
|        |    └── filtered_gene_bc_matrices
|        |                |── barcodes.tsv.gz
|        |                |── features.tsv.gz
|        |                └── matrix.mtx.gz
|        └── Mlo
|             └── filtered_gene_bc_matrices
|                         |── barcodes.tsv.gz
|                         |── features.tsv.gz
|                         └── matrix.mtx.gz
└── MZBTrajectory 
         |── Sample1
         |     └── filtered_gene_bc_matrices
         |                 |── barcodes.tsv.gz
         |                 |── features.tsv.gz
         |                 └── matrix.mtx.gz
         |── Sample2
         |     └── filtered_gene_bc_matrices
         |                 |── barcodes.tsv.gz
         |                 |── features.tsv.gz
         |                 └── matrix.mtx.gz
         └── Sample3
                └── filtered_gene_bc_matrices
                           |── barcodes.tsv.gz
                           |── features.tsv.gz
                           └── matrix.mtx.gz
