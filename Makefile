.PHONY: test-docker build-docker transpose svd

docker-name = svd-competition
pwd=$(shell pwd)

gen-data=$(shell find data/*/*200.mat)
svd-data=$(patsubst data/%.mat,data/truth/%.mat,$(gen-data))
svd:$(svd-data)

origin-data=$(shell find data/*/*.mat)
trans-data=$(patsubst data/%.mat,data/trans/%.mat,$(origin-data))
transpose:$(trans-data)

data/truth/%.mat:data/trans/%.mat
	docker run -it --rm -v $(pwd):/home/svd_competition $(docker-name) \
	python code/gen_truth.py --input_path $(dir $<) --input_name $(notdir $@) --output_path $(dir $@) --svd

data/trans/%.mat:data/%.mat
	docker run -it --rm -v $(pwd):/home/svd_competition $(docker-name) \
	python code/gen_truth.py --input_path $(dir $<) --input_name $(notdir $@) --output_path $(dir $@) --transpose

test-docker:build-docker
	docker run -it --rm -v $(pwd):/home/svd-competition $(docker-name)

build-docker:
	docker build -t $(docker-name) .