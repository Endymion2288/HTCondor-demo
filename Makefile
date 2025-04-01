# 定义关键路径和参数
DATASET_LIST   := dataset.list
MERGED_OUTPUT  := filtered_data_preCut_mu.root
TEMPLATE_FILE  := ParticleCand/src/runReadTreeX.C
ROOT_CMD       := root -l -b -q

# 从清单文件获取输入文件列表
INPUT_FILES    := $(shell cat $(DATASET_LIST))

# 生成处理后的文件路径列表
PROCESSED_FILES := $(foreach f,$(INPUT_FILES),\
					 $(dir $(f))$(basename $(notdir $(f)))-filtered/filtered_data_preCut_mu.root)
PROCESSED_SCRIPTS := $(foreach f,$(INPUT_FILES),\
					 $(dir $(f))$(basename $(notdir $(f)))-filtered/runReadTreeX.C)

# 默认目标：执行完整处理流程
all: $(MERGED_OUTPUT)

# Merged target: merge all processed files, depend on the scripts.
$(MERGED_OUTPUT): dataset_filtered.list
	@echo "Merging processed files..."
	@hadd -f $@ `cat $<`
	@echo "\n=== Merged output created: $@ ==="

# List of non-empty output files: depend on the scripts.
dataset_filtered.list: $(PROCESSED_SCRIPTS)
	@echo "Creating dataset list..."
	@echo "$(PROCESSED_SCRIPTS)"
	@rm -f $@
	@for script in $(PROCESSED_SCRIPTS); do\
		find "$$(dirname "$$script")" -regex ".*root";\
	done | grep "filtered_data_preCut_mu.root" | sort -u > $@

# Single file processing: create script and run.
# - Indicate if any product is created.
%-filtered/runReadTreeX.C: %.root
	@echo "Processing $*..."
	@mkdir -p $*-filtered/
	@sed -e "s|JOB_DATA|$<|g" $(TEMPLATE_FILE) | \
		sed -e "s|// #define RUN_JOB|#define RUN_JOB|g" > $@
	@echo "Created $@"
	@cd $*-filtered && $(ROOT_CMD) $@
	@echo "Processed $<"
	@if [ -f $*-filtered/filtered_data_preCut_mu.root ]; then \
		echo "=== $< processed into $*-filtered/filtered_data_preCut_mu.root ==="; \
	else \
		echo "=== No valid candidate found in $< ==="; \
	fi


# 清理目标
clean:
	@rm -f $(MERGED_OUTPUT)
	@for f in $(PROCESSED_FILES); do \
		dir=$$(dirname $$f); \
		rm -rf $$dir; \
		echo "Cleaned $$dir"; \
	done
	@echo "=== All processed directories cleaned ==="

.PHONY: all clean