.PHONY: test-docker build-docker

docker-name = svd-competition
pwd=$(shell pwd)

test-docker:build-docker
	docker run -it --rm -v $(pwd):/home/svd-competition $(docker-name)
build-docker:
	docker build -t $(docker-name) .