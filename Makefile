# 定义关键路径和参数
DATADIR_LIST   := datadir.list
MERGED_OUTPUT  := filtered_data_preCut_Jpsi_p6_2022_2024.root
TEMPLATE_FILE  := ParticleCand/src/runReadTreeX.C
MACRO_FILE     := ParticleCand/src/ReadTree.C
ROOT_CMD       := root -l -b -q

INPUT_DIRS        := $(shell cat $(DATADIR_LIST))
PROCESSED_SCRIPTS := $(addsuffix /runReadTreeX.C, $(INPUT_DIRS))
PROCESSED_OUTPUTS := $(addsuffix /filtered_data_preCut_mu.root, $(INPUT_DIRS))

# Phony targets: all clean
.PHONY: all clean

# Target: all
all: $(MERGED_OUTPUT)

# Target: $(MERGED_OUTPUT)
# Merge with hadd, from all the input directories
$(MERGED_OUTPUT): $(PROCESSED_OUTPUTS) $(TEMPLATE_FILE) $(DATADIR_LIST)
	@echo "Merging all the outputs..."
	@hadd -f $(MERGED_OUTPUT) $(PROCESSED_OUTPUTS)

# Target: $(PROCESSED_OUTPUTS)
# Run the script in each input directory
%/filtered_data_preCut_mu.root: %/runReadTreeX.C 
	@echo "Processing $<..."
	@cd $(dir $<) && $(ROOT_CMD) runReadTreeX.C

# Target: $(PROCESSED_SCRIPTS)
# Copy the template script to each input directory. Do substitution to create.
%/runReadTreeX.C: $(TEMPLATE_FILE) $(MACRO_FILE)
	@echo "Creating $@..."
	@sed -e 's|JOB_DATA|$*/mymultilep*.root|g' $(TEMPLATE_FILE)|\
	 sed -e 's|// #define RUN_JOB|#define RUN_JOB|g' > $@

# Target: clean
clean:
	@echo "Cleaning up..."
	@rm -f $(MERGED_OUTPUT)
	@rm -f $(PROCESSED_SCRIPTS)
	@rm -f $(PROCESSED_OUTPUTS)
	@echo "Done."