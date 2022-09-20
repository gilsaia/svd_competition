.PHONY: test-docker build-docker transpose svd simple-check simple-check-task1 simple-check-task2 simple-check-task3

docker-name = svd-competition
pwd=$(shell pwd)

simple-check:simple-check-task1 simple-check-task2 simple-check-task3

simple-check-task1:transpose
	docker run -it --rm -v $(pwd):/home/svd_competition $(docker-name) \
	python utils/check.py --simple_check --input_path data/trans/task1/ --output_path res/ --task 1

simple-check-task2:transpose
	docker run -it --rm -v $(pwd):/home/svd_competition $(docker-name) \
	python utils/check.py --simple_check --input_path data/trans/task2/ --output_path res/ --task 2

simple-check-task3:transpose
	docker run -it --rm -v $(pwd):/home/svd_competition $(docker-name) \
	python utils/check.py --simple_check --input_path data/trans/task3/ --output_path res/ --task 3

origin-data=$(shell find data/*/*.mat)
trans-data=$(patsubst data/%.mat,data/trans/%.mat,$(origin-data))
transpose:$(trans-data)

gen-data=$(shell find data/task[123]/*200.mat)
svd-data=$(patsubst data/%.mat,data/svd/%.mat,$(gen-data))
svd:$(svd-data)

data/svd/%.mat:data/trans/%.mat utils/gen_truth.py
	docker run -it --rm -v $(pwd):/home/svd_competition $(docker-name) \
	python utils/gen_truth.py --input_path $(dir $<) --input_name $(notdir $@) --output_path $(dir $@) --svd

data/trans/%.mat:data/%.mat utils/gen_truth.py
	docker run -it --rm -v $(pwd):/home/svd_competition $(docker-name) \
	python utils/gen_truth.py --input_path $(dir $<) --input_name $(notdir $@) --output_path $(dir $@) --transpose

test-docker:build-docker
	docker run -it --rm -v $(pwd):/home/svd_competition $(docker-name)

build-docker:
	docker build -t $(docker-name) .