-include .env
.PHONY: test-docker build-docker transpose svd simple-check simple-check-task1 simple-check-task2 simple-check-task3

docker-name = svd-competition
pwd=$(shell pwd)

ifeq ($(local),on)
	docker-cmd= 
else
	docker-cmd=docker run -it --rm -v $(pwd):/home/svd_competition -v $(pwd)/.wandb/.netrc:/root/.netrc $(docker-name)
endif

ifeq ($(matlab),on)
	matlab-option=--matlab
	one-matrix-cmd=matlab -nodesktop -nosplash --path code/ --path utils/ -r utils/one_matrix_test
else
	one-matrix-cmd=octave-cli --path code/ --path utils/ utils/one_matrix_test.m
endif

ifeq ($(wandb),offline)
	wandb-option=--wandb offline
endif

ifdef name
	run-name=--run_name $(name)
endif

measure:measure-task1 measure-task2 measure-task3 measure-task4

measure-task1:transpose
	$(docker-cmd) python utils/check.py --measure $(run-name) --input_path data/trans/task1/ --output_path res/ --task 1 $(matlab-option) $(wandb-option)
measure-task2:transpose
	$(docker-cmd) python utils/check.py --measure $(run-name) --input_path data/trans/task2/ --output_path res/ --task 2 $(matlab-option) $(wandb-option)
measure-task3:transpose
	$(docker-cmd) python utils/check.py --measure $(run-name) --input_path data/trans/task3/ --output_path res/ --task 3 $(matlab-option) $(wandb-option)
measure-task4:transpose
	$(docker-cmd) python utils/check.py --measure $(run-name) --input_path data/trans/task4/ --output_path res/ --task 4 $(matlab-option) $(wandb-option)

complete-check:complete-check-task1 complete-check-task2 complete-check-task3 complete-check-task4

complete-check-task1:transpose
	$(docker-cmd) python utils/check.py --complete_check --input_path data/trans/task1/ --output_path res/ --task 1 $(matlab-option)
complete-check-task2:transpose
	$(docker-cmd) python utils/check.py --complete_check --input_path data/trans/task2/ --output_path res/ --task 2 $(matlab-option)
complete-check-task3:transpose
	$(docker-cmd) python utils/check.py --complete_check --input_path data/trans/task3/ --output_path res/ --task 3 $(matlab-option)
complete-check-task4:transpose
	$(docker-cmd) python utils/check.py --complete_check --input_path data/trans/task4/ --output_path res/ --task 4 $(matlab-option)


simple-check:simple-check-task1 simple-check-task2 simple-check-task3 simple-check-task4

simple-check-task1:transpose
	$(docker-cmd) python utils/check.py --complete_check --data_len=1 --input_path data/trans/task1/ --output_path res/ --task 1 $(matlab-option)
simple-check-task2:transpose
	$(docker-cmd) python utils/check.py --complete_check --data_len=1 --input_path data/trans/task2/ --output_path res/ --task 2 $(matlab-option)
simple-check-task3:transpose
	$(docker-cmd) python utils/check.py --complete_check --data_len=1 --input_path data/trans/task3/ --output_path res/ --task 3 $(matlab-option)
simple-check-task4:transpose
	$(docker-cmd) python utils/check.py --complete_check --data_len=1 --input_path data/trans/task4/ --output_path res/ --task 4 $(matlab-option)

one-matrix-task1:
	$(docker-cmd) $(one-matrix-cmd) 1
one-matrix-task2:
	$(docker-cmd) $(one-matrix-cmd) 2
one-matrix-task3:
	$(docker-cmd) $(one-matrix-cmd) 3
one-matrix-task4:
	$(docker-cmd) $(one-matrix-cmd) 4

origin-data=$(shell find data/*/*.mat)
trans-data=$(patsubst data/%.mat,data/trans/%.mat,$(origin-data))
transpose:$(trans-data)

gen-data=$(shell find data/task[123]/*200.mat)
svd-data=$(patsubst data/%.mat,data/svd/%.mat,$(gen-data))
svd:$(svd-data)

data/svd/%.mat:data/trans/%.mat utils/gen_truth.py
	$(docker-cmd) python utils/gen_truth.py --input_path $(dir $<) --input_name $(notdir $@) --output_path $(dir $@) --svd

data/trans/%.mat:data/%.mat utils/gen_truth.py
	$(docker-cmd) python utils/gen_truth.py --input_path $(dir $<) --input_name $(notdir $@) --output_path $(dir $@) --transpose

.wandb/.netrc:
	mkdir -p .wandb
	touch .wandb/.netrc

test-docker:build-docker
	docker run -it --rm -v $(pwd):/home/svd_competition -v $(pwd)/.wandb/.netrc:/root/.netrc $(docker-name)

build-docker:.wandb/.netrc
	docker build -t $(docker-name) .