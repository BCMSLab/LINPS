#!/bin/bash

# Define variables
DATA=data
CODE=code
RESULTS=results
LOG=log

all: directories \
	$(LOG)/expression_signatures.R.Rout \
	$(LOG)/differential_expression.R.Rout \
	$(LOG)/causal_networks.R.Rout \
	$(LOG)/network_scoring.R.Rout \
	$(LOG)/similarity.R.Rout \
	$(LOG)/test_restults.out \
	$(LOG)/database.R.Rout \
	$(LOG)/test_database.out \
	$(LOG)/database_report.txt

directories: 
	mkdir -p $(DATA)
	mkdir -p $(RESULTS)
	mkdir -p $(LOG)
	@echo "Make the directories data/, results/ and log/."
	
$(LOG)/expression_signatures.R.Rout: $(CODE)/expression_signatures.R
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Download the gene expression data."
	
$(LOG)/differential_expression.R.Rout: $(CODE)/differential_expression.R \
	$(DATA)/trt_cp.rds \
	$(DATA)/trt_sh.rds
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Perform differential expression analysis."
	
$(LOG)/causal_networks.R.Rout: $(CODE)/causal_networks.R
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Download and preprocess the causal networks."
	
$(LOG)/network_scoring.R.Rout: $(CODE)/network_scoring.R \
	$(LOG)/causal_networks.R.Rout \
	$(LOG)/differential_expression.R.Rout
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Score the networks."

$(LOG)/similarity.R.Rout: $(CODE)/similarity.R \
	$(LOG)/differential_expression.R.Rout
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Calculate similarities."
	
$(LOG)/database.R.Rout: $(CODE)/database.R \
	$(LOG)/network_scoring.R.Rout
	R CMD BATCH --vanilla $< $(LOG)/$(<F).Rout
	@echo "Build the database."
	
$(LOG)/test_restults.out: tests/test-restults.R \
	$(LOG)/network_scoring.R.Rout
	Rscript -e "testthat::test_file('tests/test-restults.R')" > $(LOG)/test_restults.out 2>&1
	
$(LOG)/test_database.out: tests/test-database.R \
	$(RESULTS)/LINPS.sqlite
	Rscript -e "testthat::test_file('tests/test-database.R')" > $(LOG)/test_database.out 2>&1

$(LOG)/database_report.txt: $(CODE)/database_report.sh \
	$(RESULTS)/LINPS.sqlite
	bash $(CODE)/database_report.sh > $(LOG)/database_report.txt